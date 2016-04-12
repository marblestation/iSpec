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
import ctypes
import time
from datetime import datetime, timedelta
import numpy as np
from mpfitmodels import MPFitModel
#from continuum import fit_continuum
from abundances import write_solar_abundances, write_fixed_abundances, determine_abundances, create_free_abundances_structure
from abundances import determine_abundance_enchancements, enhance_solar_abundances
from segments import create_segments_around_lines
from atmospheres import write_atmosphere, interpolate_atmosphere_layers, valid_atmosphere_target, calculate_opacities, extrapolated_model_atmosphere_were_used
from lines import write_atomic_linelist, write_isotope_data, _get_atomic_linelist_definition, _sampling_uniform_in_velocity
from common import mkdir_p, estimate_vmic, estimate_vmac, which
from common import is_turbospectrum_support_enabled, is_spectrum_support_enabled, is_moog_support_enabled, is_synthe_support_enabled
from common import is_sme_support_enabled
from spectrum import create_spectrum_structure, convolve_spectrum, correct_velocity, resample_spectrum, read_spectrum, normalize_spectrum, create_wavelength_filter, read_spectrum, write_spectrum
from synth_sme import _sme_librrayversion, _sme_inputlinelist, _sme_inputmodel, _sme_inputabund, _sme_ionization, _sme_setvwscale, _sme_inputwaverange, _sme_opacity, _sme_transf
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
from lockfile import FileLock, LockTimeout, AlreadyLocked
from Queue import Empty
import subprocess
import shutil
import re
import glob

from scipy.signal import fftconvolve
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
    SYNTH_STEP_VRAD = 5
    SYNTH_STEP_ABUNDANCES = 0.05
    SYNTH_STEP_LOGGF = 0.01
    ###################################
    EW_STEP_TEFF = 500.
    EW_STEP_LOGG = 0.5
    EW_STEP_MH = 0.05
    EW_STEP_VMIC = 0.5


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, molecules_files=None, isotope_file=None, regions=None, waveobs_mask=None, code="spectrum", use_molecules=False, tmp_dir=None):
    """
    waveobs_mask is for SPECTRUM
    regions is for Turbospectrum
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    # Filter out lines not supported by a given synthesizer
    lcode = linelist[code+'_support'] == "T"
    linelist = linelist[lcode]
    # Limit linelist to the region asked to be synthesized
    # Provide some margin or near-by deep lines might be omitted
    margin = 2. # 2 nm
    wfilter = np.logical_and(linelist['wave_nm'] >= np.min(waveobs)-margin, linelist['wave_nm'] <= np.max(waveobs)+margin)
    linelist = linelist[wfilter]

    if code == "turbospectrum":
        return __turbospectrum_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "moog":
        return __moog_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "synthe":
        return __synthe_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "sme":
        return __sme_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose, regions=regions, timeout=timeout)
    elif code == "spectrum":
        return __spectrum_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions, waveobs_mask=waveobs_mask)


def __spectrum_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None, waveobs_mask=None):
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
    return __spectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=0.0, vsini=0.0, limb_darkening_coeff=0.00, R=0, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions)


def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, molecules_files=None, isotope_file=None, regions=None, waveobs_mask=None, code="spectrum", use_molecules=False, tmp_dir=None):
    """
    waveobs_mask is for SPECTRUM
    regions is for Turbospectrum, moog and synthe
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    # Filter out lines not supported by a given synthesizer
    lcode = linelist[code+'_support'] == "T"
    linelist = linelist[lcode]
    # Limit linelist to the region asked to be synthesized
    # Provide some margin or near-by deep lines might be omitted
    margin = 2. # 2 nm
    wfilter = np.logical_and(linelist['wave_nm'] >= np.min(waveobs)-margin, linelist['wave_nm'] <= np.max(waveobs)+margin)
    linelist = linelist[wfilter]

    if code == "turbospectrum":
        return __turbospectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "moog":
        return __moog_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "synthe":
        return __synthe_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
    elif code == "sme":
        return __sme_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose, regions=regions)
    elif code == "spectrum":
        return __spectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions, waveobs_mask=waveobs_mask)


def __spectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=0., vsini=0., limb_darkening_coeff=0., R=0, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None, waveobs_mask=None, tmp_dir=None):
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
    if len(linelist) > 500000:
        raise Exception("Linelist too big for SPECTRUM: %i (limit 500000)" % (len(linelist)))
    if fixed_abundances is None:
        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])

    if waveobs_mask is None:
        if regions is None:
            waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
            # Limit linelist
            wave_base = np.min(waveobs)
            wave_top = np.max(waveobs)
            # Provide some margin or near-by deep lines might be omitted
            margin = 2. # 2 nm
            lfilter = np.logical_and(linelist['wave_A'] >= (wave_base-margin)*10., linelist['wave_A'] <= (wave_top+margin)*10.)
            linelist = linelist[lfilter]
        else:
            waveobs_mask = _create_waveobs_mask(waveobs, regions)
            # Limit linelist
            linelist = __filter_linelist(linelist, regions)

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_fixed_abund_file = False
    remove_tmp_linelist_file = False
    remove_tmp_isotope_file = False
    if atmosphere_layers_file is None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, code="spectrum", tmp_dir=tmp_dir)
        remove_tmp_atm_file = True
    if abundances_file is None:
        abundances_file = write_solar_abundances(abundances, tmp_dir=tmp_dir)
        remove_tmp_abund_file = True
    if fixed_abundances_file is None:
        fixed_abundances_file = write_fixed_abundances(fixed_abundances, tmp_dir=tmp_dir)
        remove_tmp_fixed_abund_file = True
    if linelist_file is None:
        linelist_file = write_atomic_linelist(linelist, code="spectrum", tmp_dir=tmp_dir)
        remove_tmp_linelist_file = True
    if isotope_file is None:
        isotope_file = write_isotope_data(isotopes, tmp_dir=tmp_dir)
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

    p = Process(target=__spectrum_true_generate_spectrum, args=(process_communication_queue, waveobs, waveobs_mask, atmosphere_layers_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel), kwargs={'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose})
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


def calculate_theoretical_ew_and_depth(atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, microturbulence_vel = 2.0, atmosphere_layers_file=None, abundances_file=None, linelist_file=None, isotope_file=None, verbose=0, gui_queue=None, timeout=1800, tmp_dir=None):
    """
    """
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    supported = linelist['spectrum_support'] == "T"

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_linelist_file = False
    remove_tmp_isotope_file = False
    if atmosphere_layers_file is None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, code="spectrum", tmp_dir=tmp_dir)
        remove_tmp_atm_file = True
    if abundances_file is None:
        abundances_file = write_solar_abundances(abundances, tmp_dir=tmp_dir)
        remove_tmp_abund_file = True
    if linelist_file is None:
        linelist_file = write_atomic_linelist(linelist[supported], code="spectrum", tmp_dir=tmp_dir)
        remove_tmp_linelist_file = True
    if isotope_file is None:
        isotope_file = write_isotope_data(isotopes, tmp_dir=tmp_dir)
        remove_tmp_isotope_file = True
    nlayers = len(atmosphere_layers)
    start = np.min(linelist['wave_A'][supported]) - 0.1
    end = np.max(linelist['wave_A'][supported]) + 0.1
    num_lines = len(linelist[supported])

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
    #resulting_linelist["valid_theoretical_ew_depth"] = np.abs(output_wave - linelist['wave_A']) < 1e-5
    resulting_linelist["theoretical_ew"][supported] = np.round(output_ew, 2)
    resulting_linelist["theoretical_depth"][supported] = np.round(output_depth, 2)
    return resulting_linelist

def __enqueue_progress(process_communication_queue, v):
    process_communication_queue.put(("self.update_progress(%i)" % v))
    process_communication_queue.join()

def __spectrum_true_generate_spectrum(process_communication_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence=0., vsini=0., limb_darkening_coeff=0., R=0, nlayers=56, verbose=0):
    """
    Generate synthetic spectrum and apply macroturbulence, rotation (visini), limb darkening coeff and resolution except
    if all those parameters are set to zero, in that case the fundamental synthetic spectrum is returned.
    """
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    import synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    ## The convolution (R), rotation broadening (vsini) and macroturbulence broadening (vmac),
    ## do not seem to work as expected in the SPECTRUM code, thus we set them to zero and
    ## we use a python implementation
    #fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, nlayers, verbose, update_progress_func)
    fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel, 0, 0, 0, 0, nlayers, verbose, update_progress_func)

    segments = None
    vrad = (0,)
    fluxes = apply_post_fundamental_effects(waveobs, fluxes, segments, \
                macroturbulence=macroturbulence, vsini=vsini, \
                limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)

    process_communication_queue.put(fluxes)

def __calculate_ew_and_depth(process_communication_queue, atmosphere_model_file, linelist_file, isotope_file, abundances_file, num_lines, microturbulence_vel = 2.0, nlayers=56, start=3000, end=11000, verbose=0):
    """
    start and end in Amstrom
    """
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    import synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    output_wave, output_code, output_ew, output_depth = synthesizer.calculate_ew_and_depth(atmosphere_model_file, linelist_file, isotope_file, abundances_file, num_lines, microturbulence_vel, nlayers, start, end, verbose, update_progress_func)
    process_communication_queue.put((output_wave, output_code, output_ew, output_depth))


def apply_post_fundamental_effects(waveobs, fluxes, segments, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.20, R=500000, vrad=(0,), verbose=0):
    """
    Apply macroturbulence, rotation (visini), limb darkening coeff and resolution to already generated fundamental synthetic spectrum.
    """
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10

    spectrum = create_spectrum_structure(waveobs, fluxes)
    spectrum.sort(order=['waveobs'])

    if (macroturbulence is not None and macroturbulence > 0) or (vsini is not None and vsini > 0):
        # Build spectrum with sampling uniform in velocity (required by vmac and vsini broadening):
        wave_base = spectrum['waveobs'][0]
        wave_top = spectrum['waveobs'][-1]
        velocity_step = __determine_velocity_step(spectrum)
        waveobs_uniform_in_velocity = _sampling_uniform_in_velocity(wave_base, wave_top, velocity_step)
        fluxes_uniform_in_velocity = np.interp(waveobs_uniform_in_velocity, spectrum['waveobs'], spectrum['flux'], left=0.0, right=0.0)

        # Apply broadening
        fluxes_uniform_in_velocity = __vsini_broadening_limbdarkening(fluxes_uniform_in_velocity, velocity_step, vsini, limb_darkening_coeff)
        fluxes_uniform_in_velocity = __vmac_broadening(fluxes_uniform_in_velocity, velocity_step, macroturbulence)

        # Resample to origin wavelength grid
        fluxes = np.interp(spectrum['waveobs'], waveobs_uniform_in_velocity, fluxes_uniform_in_velocity, left=0.0, right=0.0)
        spectrum['flux'] = fluxes

    if R is not None and R > 0:
        # Convolve (here it is not needed to be with a sampling uniform in velocity, the function is capable of dealing with that)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']

    # Make sure original zeros are set to 1.0 and not modified by the previous broadening operations
    fluxes[zeros] = 1.0e-10

    if type(vrad) not in (tuple, list, np.ndarray):
        raise Exception("Velocity should be an array")

    if np.any(np.asarray(vrad) > 0):
        if len(vrad) != len(segments):
            raise Exception("Velocity should be an array with as many numbers as segments when segments are provided")
        modified = waveobs < 0 # All to false
        for velocity, segment in zip(vrad, segments):
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
            modified = np.logical_or(modified, wfilter)
            spectrum = create_spectrum_structure(waveobs[wfilter], fluxes[wfilter])
            spectrum = correct_velocity(spectrum, velocity)
            spectrum = resample_spectrum(spectrum, waveobs[wfilter], method="linear", zero_edges=True)
            fluxes[wfilter] = spectrum['flux']
        fluxes[~modified] = 1.
    return fluxes


class SynthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, isotopes, linelist_free_loggf, abundances, enhance_abundances=True, scale=None, teff=5000, logg=3.0, MH=0.0, vmic=2.0, vmac=0.0, vsini=2.0, limb_darkening_coeff=0.0, R=0, precomputed_grid_dir=None):
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
        p = [teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, R ]
        #
        self.abundances_file = None
        self.linelist_file = None
        self.isotope_file = None
        self.molecules_files = None
        self.atmosphere_layers_file = None
        super(SynthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        key = "%.0f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
        complete_key = "%.0f %.2f %.2f %.2f %.2f %.2f %.2f %i " % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()))

        # Consider new loggf
        linelist_free_loggf = self.generate_linelist_free_loggf()

        loggf_key = " ".join(map(lambda x: "%.3f" % (x), linelist_free_loggf['loggf']))
        complete_key += loggf_key
        key += loggf_key

        if len(linelist_free_loggf) > 0:
            linelist = np.hstack((self.linelist, linelist_free_loggf))
            linelist.sort(order=['wave_nm'])
        else:
            linelist = self.linelist

        # Consider new abundances as fixed
        fixed_abundances = self.free_abundances()

        abundances_key = " ".join(map(lambda x: "%.2f" % (x), fixed_abundances['Abund']))
        complete_key += abundances_key
        key += abundances_key

        ##### [start] Check precomputed (solar abundance)
        precomputed_file = str(self.precomputed_grid_dir) + "/unconvolved_steps/{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}.fits".format(int(self.teff()), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff())
        fundamental_precomputed_file = str(self.precomputed_grid_dir) + "/unconvolved_steps/{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}.fits".format(int(self.teff()), self.logg(), self.MH(), self.vmic(), 0., 0., 0.)
        if self.precomputed_grid_dir is not None and abundances_key == "" and os.path.exists(precomputed_file):
            if not self.quiet:
                print "Pre-computed:", complete_key
            precomputed = read_spectrum(precomputed_file)
            convolved_precomputed = convolve_spectrum(precomputed, self.R())

            convolved_precomputed = resample_spectrum(convolved_precomputed, self.waveobs, method="linear", zero_edges=True)
            convolved_precomputed['flux'][self.waveobs_mask == 0] = 1.
            self.last_fluxes = convolved_precomputed['flux'].copy()
            self.last_final_fluxes = convolved_precomputed['flux'].copy()

        elif self.precomputed_grid_dir is not None and abundances_key == "" and os.path.exists(fundamental_precomputed_file):
            if not self.quiet:
                print "Pre-computed (fundamental):", complete_key
            fundamental_precomputed = read_spectrum(fundamental_precomputed_file)
            fundamental_precomputed = resample_spectrum(fundamental_precomputed, self.waveobs, method="linear", zero_edges=True)
            fundamental_precomputed['flux'][self.waveobs_mask == 0] = 1.
            self.last_fluxes = fundamental_precomputed['flux']
            # Optimization to avoid too small changes in parameters or repetition
            self.cache[key] = self.last_fluxes.copy()
            self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, self.segments, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=self.R(), vrad=self.vrad(), verbose=0)
            self.last_final_fluxes[self.waveobs_mask == 0] = 1.
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
                atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH(), code=self.code)
                # Fundamental synthetic fluxes
                if self.code == "turbospectrum":
                    self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(), atmosphere_layers_file=self.atmosphere_layers_file, abundances_file=self.abundances_file, linelist_file=self.linelist_file, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, use_molecules=self.use_molecules, tmp_dir=self.tmp_dir, timeout=self.timeout)
                elif self.code == "moog":
                    self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(), atmosphere_layers_file=self.atmosphere_layers_file, abundances_file=self.abundances_file, linelist_file=self.linelist_file, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, tmp_dir=self.tmp_dir, timeout=self.timeout)
                elif self.code == "synthe":
                    self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(), atmosphere_layers_file=self.atmosphere_layers_file, abundances_file=self.abundances_file, linelist_file=self.linelist_file, molecules_files=self.molecules_files, isotope_file=self.isotope_file, regions=self.segments, verbose=0, code=self.code, tmp_dir=self.tmp_dir, timeout=self.timeout)
                elif self.code == "sme":
                    self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(), regions=self.segments, verbose=0, code=self.code, timeout=self.timeout)
                    if np.all(self.last_fluxes == 0):
                        raise Exception("SME has failed.")
                elif self.code == "spectrum":
                    self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), linelist, self.isotopes, abundances, fixed_abundances, self.vmic(),  atmosphere_layers_file=self.atmosphere_layers_file, abundances_file=self.abundances_file, linelist_file=self.linelist_file, isotope_file=self.isotope_file, waveobs_mask=self.waveobs_mask, verbose=0, tmp_dir=self.tmp_dir, timeout=self.timeout)
                    if np.all(self.last_fluxes == 0):
                        raise Exception("SPECTRUM has failed.")
                else:
                    raise Exception("Unknown code: %s" % (self.code))


                # Optimization to avoid too small changes in parameters or repetition
                self.cache[key] = self.last_fluxes.copy()

            self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, self.segments, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=self.R(), vrad=self.vrad(), verbose=0)

        return self.last_final_fluxes[self.comparing_mask]

    def fitData(self, waveobs, segments, comparing_mask, fluxes, weights=None, parinfo=None, use_errors=False, max_iterations=20, quiet=True, code="spectrum", use_molecules=False, vmic_from_empirical_relation=True, vmac_from_empirical_relation=True, tmp_dir=None, timeout=1800):
        code = code.lower()
        if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
            raise Exception("Unknown radiative transfer code: %s" % (code))

        self.timeout = timeout
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
        self.code = code
        self.use_molecules = use_molecules
        self.vmic_from_empirical_relation = vmic_from_empirical_relation
        self.vmac_from_empirical_relation = vmac_from_empirical_relation
        self.tmp_dir = tmp_dir


        # Synthesis for wavelengths with mask different from 0.0
        self.segments = segments
        self.waveobs_mask = _create_waveobs_mask(waveobs, segments)

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

        if len(self.linelist_free_loggf) == 0:
            # Only write linelist (for optimization purposes) if there is no free loggf
            if self.code == 'synthe':
                self.linelist_file, self.molecules_files = write_atomic_linelist(self.linelist, code="synthe", tmp_dir=tmp_dir)
            elif self.code == 'sme' or self.code == 'moog':
                # moog requires two files for the linelist
                # sme does not require files
                self.linelist_file = None
            else:
                # 'turbospectrum',  'spectrum
                self.linelist_file = write_atomic_linelist(self.linelist, code=self.code, tmp_dir=tmp_dir)

        if self.code == "spectrum":
            self.isotope_file = write_isotope_data(self.isotopes, tmp_dir=tmp_dir)

        # If teff, logg and MH are fixed
        if self.code != 'sme' and parinfo[0]['fixed'] and parinfo[1]['fixed'] and parinfo[2]['fixed']:
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, parinfo[0]['value'], parinfo[1]['value'], parinfo[2]['value'])
            self.atmosphere_layers_file = write_atmosphere(atmosphere_layers, parinfo[0]['value'], parinfo[1]['value'], parinfo[2]['value'], code=self.code, atmosphere_filename=None, tmp_dir=tmp_dir)

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
        if len(self.linelist_free_loggf) == 0 and self.code != 'sme' and self.code != 'moog':
            os.remove(self.linelist_file)
        if self.code == 'spectrum':
            os.remove(self.isotope_file)
        if self.code == 'synthe' and self.molecules_files is not None:
            for molecules_file in self.molecules_files:
                os.remove(molecules_file)
        # If teff, logg and MH are fixed
        if self.code != "sme" and parinfo[0]['fixed'] and parinfo[1]['fixed'] and parinfo[2]['fixed']:
            os.remove(self.atmosphere_layers_file)
        self.abundances_file = None
        self.linelist_file = None
        self.isotope_file = None
        self.molecules_files = None
        self.atmosphere_layers_file = None

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
    def vrad(self):
        vrad = []
        base = 8
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            top = base+len(self.segments)
            for i in xrange(base, top):
                vrad.append(self._parinfo[i]['value'])
        return vrad
    def free_loggf(self):
        base = 8
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            base += len(self.segments)
        loggf = []
        for i in xrange(base, base+len(self.linelist_free_loggf)):
            loggf.append(self._parinfo[i]['value'])
        return loggf
    def free_abundances(self):
        base = 8
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            base += len(self.segments)
        base += len(self.linelist_free_loggf)
        fixed_abundances = np.recarray((len(self._parinfo)-base, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])
        for i in xrange(len(self._parinfo)-base):
            fixed_abundances['code'][i] = int(self._parinfo[base+i]['parname'])
            fixed_abundances['Abund'][i] = self._parinfo[base+i]['value']
            fixed_abundances['element'][i] = ""
        return fixed_abundances

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def eMH(self): return self.m.perror[2]
    def evmic(self): return self.m.perror[3]
    def evmac(self): return self.m.perror[4]
    def evsini(self): return self.m.perror[5]
    def elimb_darkening_coeff(self): return self.m.perror[6]
    def eR(self): return self.m.perror[7]
    def evrad(self):
        base = 8
        evrad = []
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            top = base+len(self.segments)
            for i in xrange(base, top):
                evrad.append(self.m.perror[i])
        return evrad
    def efree_loggf(self):
        base = 8
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            base += len(self.segments)
        eloggf = []
        for i in xrange(base, base+len(self.linelist_free_loggf)):
            eloggf.append(self.m.perror[i])
        return eloggf
    def efree_abundances(self):
        base = 8
        if len(self._parinfo) > base and "vrad" in self._parinfo[base]['parname']:
            base += len(self.segments)
        base += len(self.linelist_free_loggf)
        eabundances = []
        for i in xrange(len(self._parinfo)-base):
            eabundances.append(self.m.perror[base+i])
        return eabundances

    def generate_linelist_free_loggf(self):
        linelist_free_loggf = self.linelist_free_loggf.copy()
        new_loggf = self.free_loggf()
        for i in xrange(len(linelist_free_loggf)):
            linelist_free_loggf['loggf'][i] = new_loggf[i]
        return linelist_free_loggf

    def transformed_free_loggf(self):
        free_loggf = {}
        free_loggf['linelist'] = self.generate_linelist_free_loggf()
        free_loggf['loggf'] = self.free_loggf()
        free_loggf['eloggf'] = self.efree_loggf()
        return free_loggf

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

        if len(self.vrad()) >= 1:
            vrad_header = "          %8s\t%8s\t%8s\t%8s" % ("wave_base","wave_top","vrad","error")
            vrad_stats = ""
            for i, (vrad, evrad, segment) in enumerate(zip(self.vrad(), self.evrad(), self.segments)):
                vrad_stats += "Segment   %8.2f\t%8.2f\t%8.2f\t%8.2f\n" % (segment['wave_base'], segment['wave_top'], vrad, evrad)
            print ""
            print vrad_header
            print vrad_stats
            print ""


        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""
        if len(transformed_abund) > 0:
            print "           ", abundances_header
            print "Abundances:", abundances_solution
            print "Ab. errors:", abundances_errors
            print ""

        transformed_free_loggf = self.transformed_free_loggf()
        transformed_linelist_free_loggf = transformed_free_loggf['linelist']
        loggf = transformed_free_loggf['loggf']
        eloggf = transformed_free_loggf['eloggf']
        if len(transformed_linelist_free_loggf) > 0:
            loggf_header = "          %8s\t%8s\t%8s\t%8s" % ("wave_base","wave_top","log(gf)","error")
            loggf_stats = ""
            for i in xrange(len(transformed_linelist_free_loggf)):
                loggf_stats += "log(gf)     %8.4f\t%8s\t%8.3f\t%8.3f\n" % (transformed_linelist_free_loggf['wave_nm'][i], transformed_linelist_free_loggf['element'][i], loggf[i], eloggf[i])
            print ""
            print loggf_header
            print loggf_stats
            print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","wchisq","rwchisq","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.wchisq, self.reduced_wchisq, self.chisq, self.reduced_chisq, self.rms)
        if extrapolated_model_atmosphere_were_used(self.modeled_layers_pack, self.teff(), self.logg(), self.MH()):
            print ""
            print "WARNING: Extrapolated model atmospheres were used for the final solution"
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status




def __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, free_abundances, linelist_free_loggf, teff_range, logg_range, MH_range, vmic_from_empirical_relation, vmac_from_empirical_relation):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 8
    free_params = [param.lower() for param in free_params]
    if "vrad" in free_params or np.any(initial_vrad != 0):
        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base+len(initial_vrad)+len(free_abundances)+len(linelist_free_loggf))]
    else:
        parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base+len(free_abundances)+len(linelist_free_loggf))]
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
    if vmic_from_empirical_relation:
        parinfo[3]['tied'] = 'estimate_vmic(p[0], p[1], p[2])'
    parinfo[3]['step'] = Constants.SYNTH_STEP_VMIC # For auto-derivatives
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [0.0, 50.0]
    #
    parinfo[4]['parname'] = "Vmac"
    parinfo[4]['value'] = initial_vmac
    parinfo[4]['fixed'] = not parinfo[4]['parname'].lower() in free_params
    if vmac_from_empirical_relation:
        parinfo[4]['tied'] = 'estimate_vmac(p[0], p[1], p[2])'
    parinfo[4]['step'] = Constants.SYNTH_STEP_VMAC # For auto-derivatives
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [0.0, 50.0]
    #
    parinfo[5]['parname'] = "Vsini"
    parinfo[5]['value'] = initial_vsini
    parinfo[5]['fixed'] = not parinfo[5]['parname'].lower() in free_params
    parinfo[5]['step'] = Constants.SYNTH_STEP_VSINI # For auto-derivatives
    parinfo[5]['limited'] = [True, True]
    parinfo[5]['limits'] = [0.0, 300.0]
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
    # VRAD
    if "vrad" in free_params or np.any(initial_vrad != 0):
        base = 8
        for i in xrange(len(initial_vrad)):
            parinfo[base+i]['parname'] = "vrad%03i" % (i)
            parinfo[base+i]['value'] = initial_vrad[i]
            parinfo[base+i]['fixed'] = not "vrad" in free_params
            parinfo[base+i]['step'] = Constants.SYNTH_STEP_VRAD # For auto-derivatives
            parinfo[base+i]['limited'] = [True, True]
            parinfo[base+i]['limits'] = [-30., 30]
    # ABUNDANCES
    if "vrad" in free_params:
        base = 8 + len(initial_vrad)
    else:
        base = 8
    for i in xrange(len(free_abundances)):
        parinfo[base+i]['parname'] = str(free_abundances['code'][i])
        parinfo[base+i]['value'] = free_abundances['Abund'][i]
        parinfo[base+i]['fixed'] = not parinfo[base+i]['parname'].lower() in free_params
        parinfo[base+i]['step'] = Constants.SYNTH_STEP_ABUNDANCES # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-30., 0.]
    # log(gf)
    if "vrad" in free_params:
        base = 8 + len(initial_vrad) + len(free_abundances)
    else:
        base = 8 + len(free_abundances)
    for i in xrange(len(linelist_free_loggf)):
        parinfo[base+i]['parname'] = str(linelist_free_loggf['wave_nm'][i])
        parinfo[base+i]['value'] = linelist_free_loggf['loggf'][i]
        parinfo[base+i]['fixed'] = not "loggf" in free_params
        parinfo[base+i]['step'] = Constants.SYNTH_STEP_LOGGF # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-100., 0.]

    return parinfo

def __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range, free_params, adjust_model_metalicity=False):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 4
    free_params = [param.lower() for param in free_params]
    #parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base)]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.]} for i in np.arange(base)]
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
    parinfo[0]['step'] = Constants.EW_STEP_TEFF # For auto-derivatives
    #parinfo[0]['mpside'] = 2
    #parinfo[0]['mpmaxstep'] = parinfo[0]['step'] * 1.5
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
    parinfo[1]['step'] = Constants.EW_STEP_LOGG # For auto-derivatives
    #parinfo[1]['mpside'] = 2
    #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    #parinfo[1]['mpmaxstep'] = parinfo[1]['step'] * 1.5
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [min_logg, max_logg]
    #
    parinfo[2]['parname'] = "Vmic"
    parinfo[2]['value'] = initial_vmic
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = Constants.EW_STEP_VMIC # For auto-derivatives
    #parinfo[2]['mpside'] = 2
    #parinfo[2]['mpmaxstep'] = parinfo[2]['step'] * 2.0
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., 50.0]
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
    parinfo[3]['parname'] = "MH"
    parinfo[3]['value'] = initial_MH
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = Constants.EW_STEP_MH # For auto-derivatives
    #parinfo[3]['mpside'] = 2
    #if not parinfo[3]['fixed']:
        #parinfo[3]['mpmaxstep'] = parinfo[3]['step'] * 1.5
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [min_MH, max_MH]

    return parinfo


def __filter_linelist(linelist, segments):
    # Provide some margin or near-by deep lines might be omitted
    margin = 2. # 2 nm
    # Build wavelength points from regions
    lfilter = None
    for region in segments:
        wave_base = region['wave_base'] - margin
        wave_top = region['wave_top'] + margin

        if lfilter is None:
            lfilter = np.logical_and(linelist['wave_A'] >= wave_base*10., linelist['wave_A'] <= wave_top*10.)
        else:
            lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave_A'] >= wave_base*10., linelist['wave_A'] <= wave_top*10.))

    if lfilter is not None:
        return linelist[lfilter]
    else:
        return linelist

# Single underscore name or it is not found from the function SynthModel
def _create_waveobs_mask(waveobs, segments):
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

def model_spectrum(spectrum, continuum_model, modeled_layers_pack, linelist, isotopes, abundances, free_abundances, linelist_free_loggf, initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=None, linemasks=None, enhance_abundances=True, scale=None, precomputed_grid_dir=None, use_errors=True, max_iterations=20, verbose=1, code="spectrum", use_molecules=False, vmic_from_empirical_relation=False, vmac_from_empirical_relation=False, tmp_dir=None, timeout=1800):
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

    * timeout is for single synthesis execution and not for the whole minimization.
    """

    if verbose or verbose == 1:
        verbose = True
        quiet = False
    else:
        verbose = False
        quiet = True

    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))


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
        # No free abundances
        free_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|S30')])
    else:
        # Add free abundances as free params
        for element in free_abundances:
            free_params.append(str(element['code']))

    atomic_dtype = _get_atomic_linelist_definition()
    if linelist_free_loggf is None:
        linelist_free_loggf = np.recarray((0, ), dtype=atomic_dtype)
    else:
        # Make sure both linelists have the same structure so that it will be possible to merge them
        linelist = linelist[[x[0] for x in atomic_dtype]]
        linelist_free_loggf = linelist_free_loggf[[x[0] for x in atomic_dtype]]
        # Add free loggf as free params
        free_params.append("loggf")

    if "vmic" in free_params and vmic_from_empirical_relation:
        vmic_from_empirical_relation = False
        logging.warn("'vmic_from_empirical_relation' changed to False because vmic is a free parameter")

    if "vmac" in free_params and vmac_from_empirical_relation:
        vmac_from_empirical_relation = False
        logging.warn("'vmac_from_empirical_relation' changed to False because vmac is a free parameter")

    if segments is None:
        segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        segments['wave_base'][0] = wave_base
        segments['wave_top'][0] = wave_top
        # Limit linelist
        # Provide some margin or near-by deep lines might be omitted
        margin = 2. # 2 nm
        lfilter = np.logical_and(linelist['wave_A'] >= (wave_base-margin)*10., linelist['wave_A'] <= (wave_top+margin)*10.)
        linelist = linelist[lfilter]
    else:
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


    teff_range = modeled_layers_pack[4]
    logg_range = modeled_layers_pack[5]
    MH_range = modeled_layers_pack[6]

    if type(initial_vrad) not in (np.ndarray, list, tuple):
        if segments is not None:
            initial_vrad = np.ones(len(segments))*initial_vrad
        else:
            initial_vrad = np.ones(1)*initial_vrad

    parinfo = __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, free_abundances, linelist_free_loggf, teff_range, logg_range, MH_range, vmic_from_empirical_relation, vmac_from_empirical_relation)

    synth_model = SynthModel(modeled_layers_pack, linelist, isotopes, linelist_free_loggf, abundances, enhance_abundances=enhance_abundances, scale=scale, precomputed_grid_dir=precomputed_grid_dir)

    #segments = None
    synth_model.fitData(waveobs, segments, comparing_mask, flux, weights=weights, parinfo=parinfo, use_errors=use_errors, max_iterations=max_iterations, quiet=quiet, code=code, use_molecules=use_molecules, vmic_from_empirical_relation=vmic_from_empirical_relation, vmac_from_empirical_relation=vmac_from_empirical_relation, tmp_dir=tmp_dir, timeout=timeout)

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
    for i, vrad in enumerate(synth_model.vrad()):
        params['vrad%04i' % (i)] = vrad

    errors = {}
    errors['teff'] = synth_model.eteff()
    errors['logg'] = synth_model.elogg()
    errors['MH'] = synth_model.eMH()
    errors['vmic'] = synth_model.evmic()
    errors['vmac'] = synth_model.evmac()
    errors['vsini'] = synth_model.evsini()
    errors['limb_darkening_coeff'] = synth_model.elimb_darkening_coeff()
    errors['R'] = synth_model.eR()
    for i, evrad in enumerate(synth_model.evrad()):
        errors['vrad%04i' % (i)] = evrad

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
    for i, evrad in enumerate(synth_model.evrad()):
        errors['vrad%04i' % (i)] *= error_scale_factor
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

    free_loggf = synth_model.transformed_free_loggf()

    return spectrum, synth_spectrum, params, errors, free_abundances, free_loggf, status, stats_linemasks



class EquivalentWidthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, abundances, teff=5000, logg=3.0, MH=0.0, vmic=2.0, adjust_model_metalicity=False, enhance_abundances=True, scale=None):
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

        self.min_MH = np.min(modeled_layers_pack[6])
        self.max_MH = np.max(modeled_layers_pack[6])
        #
        super(EquivalentWidthModel, self).__init__(p)



    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        key = "%.0f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
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
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH(), code=self.code)
            if self.fe1_filter is None or self.fe2_filter is None:
                ignore = np.ones(len(self.linemasks)) # Do not ignore any line since it's the first execution and it has not been done any selection
            else:
                ignore = np.zeros(len(self.linemasks))
                ignore[np.where(np.logical_or(self.fe1_filter, self.fe2_filter))[0]] = 1.0 # Do not ignore selected fe1/2 lines

            spec_abund, absolute_abund, x_over_h, x_over_fe = determine_abundances(atmosphere_layers, \
                    self.teff(), self.logg(), self.MH(), self.linemasks, self.abundances, microturbulence_vel = self.vmic(), \
                    ignore=ignore, verbose=0, code=self.code, tmp_dir=self.tmp_dir, \
                    enhance_abundances=self.enhance_abundances, scale=self.scale)


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
            #self.outliers_detection = None # Don't identify and filter outliers
            self.select_good_lines(x_over_h)

        values_to_evaluate = []
        fitted_lines_params = []
        selected_x_over_h = []

        import statsmodels.api as sm
        ### Temperature
        ## y = mx + c
        #x = self.linemasks['lower_state_eV'][self.fe1_filter]
        #y = x_over_h[self.fe1_filter]
        x = self.linemasks['lower_state_eV'][np.logical_or(self.fe1_filter, self.fe2_filter)]
        y = x_over_h[np.logical_or(self.fe1_filter, self.fe2_filter)]
        unknown = np.isnan(y)
        x = x[~unknown]
        y = y[~unknown]
        if len(x) < 2:
            raise Exception("Not enough abundances were calculated")
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m1 = linear_model.params[0]
        self.c1 = linear_model.params[1]
        self.fe1 = np.nanmedian(x_over_h[self.fe1_filter])
        self.fe1_std = np.nanstd(x_over_h[self.fe1_filter])
        ##self.fe1 = np.median(linear_model.fittedvalues)
        ##self.fe1_std = np.std(linear_model.fittedvalues)
        #print "Fe 1", np.median(linear_model.fittedvalues), np.nanmedian(x_over_h[self.fe1_filter])
        #print "    ", np.std(linear_model.fittedvalues), np.nanstd(x_over_h[self.fe1_filter])
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.scatter(x, y)
        #plt.plot(x, self.m1*x + self.c1)
        #plt.xlabel("Lower state (eV)")
        #plt.ylabel("[Fe/H]")
        #plt.grid()
        #plt.show()

        ### Vmic
        ## y = mx + c
        #x = self.linemasks['ewr'][self.fe1_filter]
        #y = x_over_h[self.fe1_filter]
        x = self.linemasks['ewr'][np.logical_or(self.fe1_filter, self.fe2_filter)]
        y = x_over_h[np.logical_or(self.fe1_filter, self.fe2_filter)]
        unknown = np.isnan(y)
        x = x[~unknown]
        y = y[~unknown]
        if len(x) < 2:
            raise Exception("Not enough abundances were calculated")
        x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
        linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
        self.m2 = linear_model.params[0]
        self.c2 = linear_model.params[1]
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.scatter(x, y)
        #plt.plot(x, self.m2*x + self.c2)
        #plt.xlabel("Reduced EW")
        #plt.ylabel("[Fe/H]")
        #plt.grid()
        #plt.show()

        ### Fe2
        ## y = mx + c
        x = self.linemasks['ewr'][self.fe2_filter]
        if len(x) > 1:
            y = x_over_h[self.fe2_filter]
            unknown = np.isnan(y)
            x = x[~unknown]
            y = y[~unknown]
            if len(x) < 2:
                raise Exception("Not enough abundances were calculated")
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
            self.fe2 = np.nanmedian(x_over_h[self.fe2_filter])
            self.fe2_std = np.nanstd(x_over_h[self.fe2_filter])
            ##self.fe2 = np.median(linear_model.fittedvalues)
            ##self.fe2_std = np.std(linear_model.fittedvalues)
            #print "Fe 2", np.median(linear_model.fittedvalues), np.nanmedian(x_over_h[self.fe2_filter])
            #print "    ", np.std(linear_model.fittedvalues), np.nanstd(x_over_h[self.fe2_filter])
            #import matplotlib.pyplot as plt
            #plt.scatter(x, y)
            #plt.plot(x, m2*x + c2)
            #plt.show()
        else:
            self.fe2 = np.nanmedian(x_over_h[self.fe2_filter])
            self.fe2_std = np.nanstd(x_over_h[self.fe2_filter])

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
        #self.select_good_lines(x_over_h) # Modifies self.lines_for_teff and self.lines_for_vmic

        return 0


    def select_good_lines(self, x_over_h, strict=True):
        """
            Modifies self.fe1_filter and self.fe2_filter
        """
        # Out of range
        unknown = np.isnan(x_over_h)
        bad = np.logical_or(x_over_h > 1.0, x_over_h < -5)
        bad = np.logical_or(bad, unknown)
        #### Line selection
        fe1_filter = self.linemasks['element'] == "Fe 1"
        fe2_filter = self.linemasks['element'] == "Fe 2"

        if self.outliers_detection not in ['robust', 'sigma_clipping']:
            strict = False
        else:
            scrict = True

        if strict and len(np.where(~bad)[0]) > 1:
            # Outliers
            import statsmodels.api as sm
            # Do not use NaN values but allow the use of out of range values since
            # they could come back to normal values later on
            x = self.linemasks['lower_state_eV'][~unknown]
            y = x_over_h[~unknown]
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            if self.outliers_detection == "robust":
                # RLM (Robust least squares)
                # Huber's T norm with the (default) median absolute deviation scaling
                # - http://en.wikipedia.org/wiki/Huber_loss_function
                # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
                huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
                linear_model = huber_t.fit()
                reject_filter1 = linear_model.weights < self.outliers_weight_limit
                #reject_filter1 = np.logical_or(reject_filter1, bad) # Done later
            elif self.outliers_detection == "sigma_clipping":
                linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
                m1 = linear_model.params[0]
                c1 = linear_model.params[1]
                corrected_y = y - (m1*x + c1)
                sigma = np.std(corrected_y)
                reject_filter1 = np.logical_or(corrected_y > + self.sigma_level*sigma, corrected_y < -self.sigma_level*sigma)
            else:
                logging.warn("Unknown outlier detection technique: %s" % (self.outliers_detection))
                reject_filter1 = np.asarray([False]*len(x))
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['lower_state_eV'], x_over_h)
            #plt.scatter(self.linemasks['lower_state_eV'][reject_filter1], x_over_h[reject_filter1], color="red")
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
            # Do not use NaN values but allow the use of out of range values since
            # they could come back to normal values later on
            x = self.linemasks['ewr'][~unknown]
            y = x_over_h[~unknown]
            x_c = sm.add_constant(x, prepend=False) # Add a constant (1.0) to have a parameter base
            if self.outliers_detection == "robust":
                # RLM (Robust least squares)
                # Huber's T norm with the (default) median absolute deviation scaling
                # - http://en.wikipedia.org/wiki/Huber_loss_function
                # - options are LeastSquares, HuberT, RamsayE, AndrewWave, TrimmedMean, Hampel, and TukeyBiweight
                huber_t = sm.RLM(y, x_c, M=sm.robust.norms.HuberT())
                linear_model = huber_t.fit()
                reject_filter2 = linear_model.weights < self.outliers_weight_limit
                #reject_filter2 = np.logical_or(reject_filter2, bad) # Done later
            elif self.outliers_detection == "sigma_clipping":
                linear_model = sm.OLS(y, x_c).fit() # Ordinary Least Square
                m2 = linear_model.params[0]
                c2 = linear_model.params[1]
                corrected_y = y - (m1*x + c1)
                sigma = np.std(corrected_y)
                reject_filter2 = np.logical_or(corrected_y > self.sigma_level*sigma, corrected_y < -self.sigma_level*sigma)
            else:
                logging.warn("Unknown outlier detection technique: %s" % (self.outliers_detection))
                reject_filter2 = np.asarray([False]*len(x))
            #import matplotlib.pyplot as plt
            #plt.scatter(self.linemasks['ewr'], x_over_h)
            #plt.scatter(self.linemasks['ewr'][reject_filter2], x_over_h[reject_filter2], color="red")
            #plt.show()

            reject_filter_tmp = np.logical_or(reject_filter1, reject_filter2)
            known_idx = np.where(~unknown)[0]

            # unknown abundances where excluded, recover them and keep the array size
            # coherent
            reject_filter = unknown.copy()
            reject_filter[known_idx] = reject_filter_tmp

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



    def fitData(self, linemasks, outliers_detection='robust', sigma_level=3, outliers_weight_limit=0.90, parinfo=None, max_iterations=20, quiet=True, code="spectrum", tmp_dir=None):
        base = 3
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        code = code.lower()
        if code not in ['spectrum', 'turbospectrum', 'moog', 'width']:
            raise Exception("Unknown radiative transfer code: %s" % (code))


        self.code = code
        self.tmp_dir = tmp_dir
        self.outliers_detection = outliers_detection
        self.sigma_level = sigma_level
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
        if code == "moog":
            chisq_limit = None
        else:
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
        if extrapolated_model_atmosphere_were_used(self.modeled_layers_pack, self.teff(), self.logg(), MH):
            print ""
            print "WARNING: Extrapolated model atmospheres were used for the final solution"
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status


def model_spectrum_from_ew(linemasks, modeled_layers_pack, abundances, initial_teff, initial_logg, initial_MH, initial_vmic, free_params=["teff", "logg", "vmic"], adjust_model_metalicity=False, enhance_abundances=True, scale=None, max_iterations=20, outliers_detection='robust', sigma_level=3, outliers_weight_limit=0.90, code="spectrum", tmp_dir=None):
    """
    - outlier_detection:
        - 'robust': Fit a robust least square linear model, outliers_weight_limit will be use as a threshold. If it is set to zero, no outliers are filtered.
        - 'sigma_clipping': Fit a tradition least square linear model and filter X times the standard deviation (sigma_level)
    - If enhance_abundances is True, alpha elements and CNO abundances will be scaled
      depending on the metallicity.
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'width']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    teff_range = modeled_layers_pack[4]
    logg_range = modeled_layers_pack[5]
    MH_range = modeled_layers_pack[6]

    # Do not allow users to set free MH in free_params to avoid confusions
    # because metallicity is always free in this method, what we make by including MH in free_params
    # is turning on the adjustment in the metallicity models
    if "MH" in free_params or "mh" in free_params:
        raise Exception("Metallicity cannot be a free parameter!")

    if adjust_model_metalicity:
        free_params.append("MH")

    parinfo = __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range, free_params, adjust_model_metalicity=adjust_model_metalicity)


    EW_model = EquivalentWidthModel(modeled_layers_pack, abundances, MH=initial_MH, adjust_model_metalicity=adjust_model_metalicity, \
                                        enhance_abundances=enhance_abundances, scale=scale)

    lfilter = linemasks['element'] == "Fe 1"
    lfilter = np.logical_or(lfilter, linemasks['element'] == "Fe 2")
    linemasks = linemasks[lfilter]
    EW_model.fitData(linemasks, parinfo=parinfo, max_iterations=max_iterations, quiet=False, outliers_detection=outliers_detection, sigma_level=sigma_level, outliers_weight_limit=outliers_weight_limit, code=code, tmp_dir=tmp_dir)
    print "\n"
    EW_model.print_solution()

    # Calculate MH
    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values
    MH = EW_model.fe1
    eMH = EW_model.fe1_std

    used_linemasks = EW_model.linemasks[np.logical_or(EW_model.fe1_filter, EW_model.fe2_filter)]

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

    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks




def __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, code="spectrum", use_molecules=False, tmp_dir=None, locked=False):
    multiprocessing.current_process().daemon=False


    if valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        if not locked:
            lock = FileLock(filename_out+".lock")
            try:
                lock.acquire(timeout=-1)    # Don't wait
            except (LockTimeout, AlreadyLocked) as e:
                # Some other process is computing this spectrum, do not continue
                print "Skipping", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, "already locked"
                return None

        try:
            print "[started]", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution
            # Prepare atmosphere model
            atmosphere_layers = interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH, code=code)
            # Synthesis
            synth_spectrum = create_spectrum_structure(wavelengths)
            synth_spectrum['flux'] = generate_spectrum(synth_spectrum['waveobs'], \
                    atmosphere_layers, teff, logg, MH, atomic_linelist, isotopes, solar_abundances, \
                    fixed_abundances=None, microturbulence_vel = vmic, \
                    macroturbulence=vmac, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
                    R=resolution, regions=segments, verbose=0, \
                    code=code, use_molecules=use_molecules, tmp_dir=tmp_dir)
            # FITS
            write_spectrum(synth_spectrum, filename_out)
            print "[finished]", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution
        finally:
            if not locked: # Not locked in this function
                lock.release()
    else:
        raise Exception("Not valid: %i %.2f %.2f" % (teff, logg, MH))



def precompute_synthetic_grid(output_dirname, ranges, wavelengths, to_resolution, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, enhance_abundances=True, scale=None, segments=None, number_of_processes=1, code="spectrum", use_molecules=False, tmp_dir=None):
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
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

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
        vsini = 0.0 # This can be modified after synthesis if needed
        limb_darkening_coeff = 0.00 # This can be modified after synthesis if needed
        resolution = 0 # This can be modified after synthesis if needed
        # For each reference point, calculate also the variations that iSpec will perform in the first iteration
        steps =   ( # Final unconvolved spectra where vmic/vmac are free and do not follow vmic/vmac empirical relations
                    (teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff+Constants.SYNTH_STEP_TEFF, logg, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg+Constants.SYNTH_STEP_LOGG, MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH+Constants.SYNTH_STEP_MH, vmic, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic+Constants.SYNTH_STEP_VMIC, vmac, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac+Constants.SYNTH_STEP_VMAC, vsini, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac, vsini+Constants.SYNTH_STEP_VSINI, limb_darkening_coeff),
                    (teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff+Constants.SYNTH_STEP_LIMB_DARKENING_COEFF),
                    # Final unconvolved spectra where vmic/vmac are not free and do follow vmic/vmac empirical relations
                    (teff+Constants.SYNTH_STEP_TEFF, logg, MH, estimate_vmic(teff+Constants.SYNTH_STEP_TEFF, logg, MH), estimate_vmac(teff+Constants.SYNTH_STEP_TEFF, logg, MH), vsini, limb_darkening_coeff),
                    (teff, logg+Constants.SYNTH_STEP_LOGG, MH, estimate_vmic(teff, logg+Constants.SYNTH_STEP_LOGG, MH), estimate_vmac(teff, logg+Constants.SYNTH_STEP_LOGG, MH), vsini, limb_darkening_coeff),
                    (teff, logg, MH+Constants.SYNTH_STEP_MH, estimate_vmic(teff, logg, MH+Constants.SYNTH_STEP_MH), estimate_vmac(teff, logg, MH+Constants.SYNTH_STEP_MH), vsini, limb_darkening_coeff),
                    # Fundamental spectra when vmic is free
                    (teff, logg, MH, vmic, 0., 0., 0.),
                    (teff+Constants.SYNTH_STEP_TEFF, logg, MH, vmic, 0., 0., 0.),
                    (teff, logg+Constants.SYNTH_STEP_LOGG, MH, vmic, 0., 0., 0.),
                    (teff, logg, MH+Constants.SYNTH_STEP_MH, vmic, 0., 0., 0.),
                    (teff, logg, MH, vmic+Constants.SYNTH_STEP_VMIC, 0., 0., 0.),
                    # Fundamental spectra when vmic is fixed and follows vmic empirical relation
                    (teff, logg, MH, estimate_vmic(teff, logg, MH), 0., 0., 0.),
                    (teff+Constants.SYNTH_STEP_TEFF, logg, MH, estimate_vmic(teff+Constants.SYNTH_STEP_TEFF, logg, MH), 0., 0., 0.),
                    (teff, logg+Constants.SYNTH_STEP_LOGG, MH, estimate_vmic(teff, logg+Constants.SYNTH_STEP_LOGG, MH), 0., 0., 0.),
                    (teff, logg, MH+Constants.SYNTH_STEP_MH, estimate_vmic(teff, logg, MH+Constants.SYNTH_STEP_MH), 0., 0., 0.))

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

                lock = FileLock(filename_out+".lock")
                try:
                    lock.acquire(timeout=-1)    # Don't wait
                except (LockTimeout, AlreadyLocked) as e:
                    # Some other process is computing this spectrum, do not continue
                    print "Skipping", teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, "already locked"
                    continue

                try:
                    tcheck = default_timer()
                    # Validate parameters
                    __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, abundances, code=code, use_molecules=use_molecules, tmp_dir=tmp_dir, locked=True)
                    elapsed = default_timer() - tcheck

                    print "-----------------------------------------------------"
                    print "Remaining time:"
                    print "\t", (num_spec-i)*elapsed, "seconds"
                    print "\t", (num_spec-i)*(elapsed/60), "minutes"
                    print "\t", (num_spec-i)*(elapsed/(60*60)), "hours"
                    print "\t", (num_spec-i)*(elapsed/(60*60*24)), "days"
                    print "-----------------------------------------------------"
                finally:
                    lock.release()

            else:
                pool.apply_async(__generate_synthetic_fits, [filename_out, wavelengths, segments, teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, resolution, modeled_layers_pack, atomic_linelist, isotopes, abundances], kwds={'code': code, 'use_molecules': use_molecules, 'tmp_dir':tmp_dir, 'locked':False})
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
        # Only use the first spectra generated for each combination
        vmic = estimate_vmic(teff, logg, MH)
        vmac = estimate_vmac(teff, logg, MH)
        vsini = 0.0
        limb_darkening_coeff = 0.00
        resolution = 0
        reference_filename_out = "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}".format(int(teff), logg, MH, vmic, vmac, vsini, limb_darkening_coeff) + ".fits"
        reference_list.add_row((reference_filename_out, int(teff), logg, MH, vmic, vmac, vsini, limb_darkening_coeff))

        # Spectra in the grid is convolved to the specified resolution for fast comparison
        print "Quick grid:", reference_filename_out
        spectrum = read_spectrum(fits_dir + reference_filename_out)
        convolved_spectrum = convolve_spectrum(spectrum, to_resolution)

        if reference_grid is None:
            reference_grid = convolved_spectrum['flux']
        else:
            reference_grid = np.vstack((reference_grid, convolved_spectrum['flux']))

    ascii.write(reference_list, reference_list_filename, delimiter='\t')
    # Generate FITS file with grid for fast comparison
    primary_hdu = fits.PrimaryHDU(reference_grid)
    wavelengths_hdu = fits.ImageHDU(wavelengths, name="WAVELENGTHS")
    params_bintable_hdu = fits.BinTableHDU(reference_list.as_array(), name="PARAMS")
    fits_format = fits.HDUList([primary_hdu, wavelengths_hdu, params_bintable_hdu])
    fits_format.writeto(reference_grid_filename, clobber=True)



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
    initial_vsini = 0.0
    initial_limb_darkening_coeff = 0.00

    reference_grid_filename = precomputed_dir + "/reference_grid_%i.fits" % resolution
    if not os.path.exists(reference_grid_filename):
        logging.warn("Pre-computed grid does not exists for R = %i" % resolution)
    else:
        try:
            grid = fits.open(reference_grid_filename)
            grid_waveobs = np.asarray(grid['WAVELENGTHS'].data, dtype=float)
            resampled_spectrum = resample_spectrum(spectrum, grid_waveobs, method="linear")

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


def __turbospectrum_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, use_molecules=False, tmp_dir=None, timeout=1800):
    return __turbospectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)

def __turbospectrum_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, use_molecules=False, tmp_dir=None, timeout=1800):
    if not is_turbospectrum_support_enabled():
        raise Exception("Turbospectrum support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    turbospectrum_dir = ispec_dir + "/synthesizer/turbospectrum/"
    turbospectrum_data = turbospectrum_dir + "/DATA/"
    turbospectrum_bsyn_lu = turbospectrum_dir + "bin/bsyn_lu"
    molecules_dir = ispec_dir + "input/linelists/turbospectrum/molecules/"

    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top
    else:
        global_wave_base = np.min(regions['wave_base'])
        global_wave_top = np.max(regions['wave_top'])

    # TODO: decide how to stablish wave_step
    waveobs = waveobs.copy()
    waveobs.sort()
    #wave_step = waveobs[1] - waveobs[0] # 0.001
    wave_step = np.max((0.001, np.min(waveobs[1:] - waveobs[:-1])))

    # Limit linelist
    linelist = __filter_linelist(linelist, regions)


    # Turbospectrum is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    atom_abundances = abundances[abundances['code'] <= 92]
    if len(atom_abundances) != 92:
        raise Exception("No abundances for all 92 elements!")
    efilter = np.logical_and(atom_abundances['code'] != 1, atom_abundances['code'] != 2)
    atom_abundances['Abund'][efilter] += MH

    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        atom_abundances = atom_abundances.copy()
        for fixed_abundance in fixed_abundances:
            index = np.where(atom_abundances['code'] == fixed_abundance['code'])[0]
            atom_abundances['Abund'][index] = fixed_abundance['Abund']

    radius = atmosphere_layers[0][-1]
    nvalues = len(atmosphere_layers[0])
    if nvalues == 11 and radius > 2.0: # Compare to 2.0 instead of 1.0 to avoid floating point imprecisions
        logging.info("Spherical model atmosphere with radius %.2e cm" % (radius))
        spherical_model = True
    else:
        spherical_model = False

    # Turbospectrum cannot compute in a single run a big chunk of wavelength so
    # we split the computation in several pieces
    max_segment = 100. # nm
    if (global_wave_top - global_wave_base)/wave_step > max_segment/wave_step:
        segment_wave_base = np.arange(global_wave_base, global_wave_top, max_segment)
        segments = np.recarray((len(segment_wave_base),),  dtype=[('wave_base', float), ('wave_top', float)])
        segments['wave_base'] = segment_wave_base
        segments['wave_top'] = segment_wave_base + max_segment - wave_step
        segments['wave_top'][-1] = global_wave_top # Last segment should not over pass the original global limits
    else:
        segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        segments['wave_base'][0] = global_wave_base
        segments['wave_top'][0] = global_wave_top


    remove_tmp_atm_file = False
    remove_tmp_linelist_file = False
    if atmosphere_layers_file is None:
        remove_tmp_atm_file = True
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=atmosphere_layers_file, code="turbospectrum", tmp_dir=tmp_dir)

    opacities_file = calculate_opacities(atmosphere_layers_file, atom_abundances, MH, microturbulence_vel, global_wave_base-10, global_wave_top+10, wave_step, verbose=verbose, opacities_filename=None, tmp_dir=tmp_dir)

    if linelist_file is None:
        remove_tmp_linelist_file = True
        linelist_filename = write_atomic_linelist(linelist, linelist_filename=linelist_file, code="turbospectrum", tmp_dir=tmp_dir)
    else:
        linelist_filename = linelist_file


    synth_fluxes = []
    synth_waveobs = []
    non_positive_result = False
    for segment in segments:
        wave_base = segment['wave_base']
        wave_top = segment['wave_top']
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
        out.close()
        synth_spectrum_filename = out.name

        # Temporary dir
        tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
        os.symlink(turbospectrum_data, tmp_execution_dir+"/DATA")
        os.symlink(molecules_dir, tmp_execution_dir+"/molecules")
        previous_cwd = os.getcwd()
        os.chdir(tmp_execution_dir)

        command = turbospectrum_bsyn_lu
        command_input = "'LAMBDA_MIN:'  '"+str(wave_base*10.)+"'\n"
        command_input += "'LAMBDA_MAX:'  '"+str(wave_top*10.)+"'\n"
        command_input += "'LAMBDA_STEP:' '"+str(wave_step*10.)+"'\n"
        for region in regions:
            command_input += "'LAMBDA_MIN:'  '"+str(region['wave_base']*10.)+"'\n"
            command_input += "'LAMBDA_MAX:'  '"+str(region['wave_top']*10.)+"'\n"
        command_input += "'INTENSITY/FLUX:' 'Flux'\n"
        command_input += "'COS(THETA)    :' '1.00'\n"
        command_input += "'ABFIND        :' '.false.'\n"
        command_input += "'MODELOPAC:' '"+opacities_file+"'\n"
        command_input += "'RESULTFILE :' '"+synth_spectrum_filename+"'\n"
        #command_input += "'METALLICITY:'    '"+str(MH)+"'\n"
        command_input += "'METALLICITY:'    '0.00'\n" # We have done the abundance changes already
        command_input += "'ALPHA/Fe   :'    '0.00'\n"
        command_input += "'HELIUM     :'    '0.00'\n"
        command_input += "'R-PROCESS  :'    '0.00'\n"
        command_input += "'S-PROCESS  :'    '0.00'\n"
        #command_input += "'INDIVIDUAL ABUNDANCES:'   '1'\n"
        #command_input += "3  1.05\n"
        command_input += "'INDIVIDUAL ABUNDANCES:'   '"+str(len(atom_abundances))+"'\n"
        for atom_abundance in atom_abundances:
            abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
            command_input +=  "%i  %.2f\n" % (atom_abundance['code'], abund)
        #command_input += "'ISOTOPES : ' '2'\n"
        #command_input += "3.006  0.075\n"
        #command_input += "3.007  0.925\n"
        command_input += "'ISOTOPES : ' '"+str(len(isotopes))+"'\n"
        for isotope in isotopes:
            command_input += "%i.%03i  %.3f\n" % (isotope['atomic_code'], isotope['mass_number'], isotope['relative_abundance_in_the_solar_system'])

        num_molecules_files = 0
        molecules = ""
        if use_molecules:
            for filename in glob.glob("molecules/*.bsyn"):
                name, file_wave_base, file_wave_top = re.match("(.*)_(\d+)-(\d+)\.bsyn", os.path.basename(filename)).groups()
                file_wave_base = float(file_wave_base)
                file_wave_top = float(file_wave_top)
                if (file_wave_base >= wave_base and file_wave_top <= wave_top) or \
                        (wave_base >= file_wave_base and wave_base <= file_wave_top ) or \
                        (wave_top >= file_wave_base and wave_top <= file_wave_top ):
                    molecules += filename + "\n"
                    num_molecules_files += 1
            command_input += "'NFILES   :' '%i'\n" % (2 + num_molecules_files)
            command_input += molecules
        else:
            command_input += "'NFILES   :' '2'\n"
        command_input += "DATA/Hlinedata\n"
        command_input += linelist_filename + "\n"
        if spherical_model:
            command_input += "'SPHERICAL:'  'T'\n"
        else:
            command_input += "'SPHERICAL:'  'F'\n"
        command_input += "  30\n"
        command_input += "  300.00\n"
        command_input += "  15\n"
        command_input += "  1.30\n"

        # If timeout command exists in PATH, then use it to control turbospectrum execution time (it might get blocked sometimes)
        if which("timeout") is not None:
            command = "timeout %i " % (timeout) + command

        if verbose == 1:
            proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
        else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)

        # wait for the process to terminate
        out, err = proc.communicate(input=command_input)
        errcode = proc.returncode

        if errcode == 124: # TIMEOUT
            logging.error("A timeout has occurred in the turbospectrum synthesis process.")
            raise Exception("Timeout: Synthesis failed!")

        os.chdir(previous_cwd)

        try:
            data = np.loadtxt(synth_spectrum_filename)
            if len(data) == 0:
                raise Exception()
        except:
            print out
            sys.stdout.flush()
            raise Exception("Synthesis failed!")
        #synth_waveobs_tmp = np.linspace(wave_base, wave_top, len(synth_fluxes_tmp)) # Not exactly identical to turbospectrum wavelengths
        synth_waveobs_tmp = data[:,0] / 10. # Armstrong to nm
        synth_waveobs = np.hstack((synth_waveobs, synth_waveobs_tmp))

        synth_fluxes_tmp = data[:,1]
        if len(synth_fluxes_tmp) > 1 and np.isnan(synth_fluxes_tmp[-1]):
            synth_fluxes_tmp[-1] = synth_fluxes_tmp[-2] # Turbospectrum bug with gfortran, last flux is always NaN
        synth_fluxes = np.hstack((synth_fluxes, synth_fluxes_tmp))

        os.remove(synth_spectrum_filename)
        shutil.rmtree(tmp_execution_dir)

    os.remove(opacities_file)
    if remove_tmp_atm_file:
        os.remove(atmosphere_layers_file)
    if remove_tmp_linelist_file:
        os.remove(linelist_filename)

    segments = None
    vrad = (0,)
    synth_fluxes = apply_post_fundamental_effects(synth_waveobs, synth_fluxes, segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)

    synth_spectrum = create_spectrum_structure(synth_waveobs, synth_fluxes)
    synth_spectrum.sort(order=['waveobs'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(synth_spectrum['waveobs'], waveobs):
        synth_spectrum = resample_spectrum(synth_spectrum, waveobs, method="linear", zero_edges=True)

    return synth_spectrum['flux']


def __moog_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, regions=None, tmp_dir=None, timeout=1800):
    return __moog_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, tmp_dir=tmp_dir, timeout=timeout)

def __moog_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, tmp_dir=None, timeout=1800):
    if not is_moog_support_enabled():
        raise Exception("MOOG support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    moog_dir = ispec_dir + "/synthesizer/moog/"
    moog_executable = moog_dir + "MOOGSILENT"

    #molecules = linelist['molecule'] == "T"
    #if len(np.where(molecules)[0]) > 0:
        #logging.warn("The molecules have been removed from the linelist because if not MOOG is unreasonably slow")
        #linelist = linelist[~molecules]


    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top
    else:
        global_wave_base = np.min(regions['wave_base'])
        global_wave_top = np.max(regions['wave_top'])

    # TODO: decide how to stablish wave_step
    waveobs = waveobs.copy()
    waveobs.sort()
    #wave_step = waveobs[1] - waveobs[0] # 0.001
    wave_step = np.max((0.001, np.min(waveobs[1:] - waveobs[:-1])))

    # Limit linelist
    linelist = __filter_linelist(linelist, regions)

    # MOOG does not have hard-coded lines, we use a external file:
    #hydrogen_lines_file = moog_dir + "/hydrogen_moog_lines.txt"
    #hydrogen_lines = ascii.read(hydrogen_lines_file, names=["wave_A", "spectrum_moog_species", "lower_state_eV", "loggf"])
    hydrogen_lines_file = moog_dir + "DATA/Hlinedata"
    hydrogen_lines = ascii.read(hydrogen_lines_file, names=["wave_A", "spectrum_moog_species", "lower_state_eV", "loggf", "designation"])

    # MOOG is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    abundances = abundances.copy()
    efilter = np.logical_and(abundances['code'] != 1, abundances['code'] != 2)
    efilter = np.logical_and(efilter, abundances['code'] <= 92)
    abundances['Abund'][efilter] += MH

    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        for fixed_abundance in fixed_abundances:
            index = np.where(abundances['code'] == fixed_abundance['code'])[0]
            abundances['Abund'][index] = fixed_abundance['Abund']

    synth_fluxes = []
    synth_waveobs = []
    for i, region in enumerate(regions):
        # It is better not to synthesize in a single run a big chunk of wavelength so
        # we split the computation in several pieces
        max_segment = 100. # nm
        if (region['wave_top'] - region['wave_base'])/wave_step > max_segment/wave_step:
            segment_wave_base = np.arange(region['wave_base'], region['wave_top'], max_segment)
            segments = np.recarray((len(segment_wave_base),),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'] = segment_wave_base
            segments['wave_top'] = segment_wave_base + max_segment - wave_step
            segments['wave_top'][-1] = region['wave_top'] # Last segment should not over pass the original global limits
        else:
            segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'][0] = region['wave_base']
            segments['wave_top'][0] = region['wave_top']

        for segment in segments:
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']

            tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
            os.makedirs(tmp_execution_dir+"/DATA/")
            atmosphere_filename = tmp_execution_dir + "/model.in"
            linelist_filename = tmp_execution_dir + "/lines.in"
            barklem_linelist_filename = tmp_execution_dir + "/DATA/Barklem.dat"
            stronglinelist_file = tmp_execution_dir + "/stronglines.in"

            if atmosphere_layers_file is None:
                atmosphere_filename = write_atmosphere(atmosphere_layers, teff, logg, MH, code="moog", atmosphere_filename=atmosphere_filename, tmp_dir=tmp_dir)
            else:
                shutil.copyfile(atmosphere_layers_file, atmosphere_filename)

            # Write (MOOG requires a separate file for damping coeff if we want to provide rad coeff and alpha from ABO theory)
            wfilter = np.logical_and(linelist['wave_nm'] >= wave_base, linelist['wave_nm'] <= wave_top)
            linelist_filename = write_atomic_linelist(linelist[wfilter], linelist_filename=linelist_filename, code="moog", tmp_dir=tmp_dir)
            barklem_linelist_filename = write_atomic_linelist(linelist[wfilter], linelist_filename=barklem_linelist_filename, code="moog_barklem", tmp_dir=tmp_dir)

            molecules = linelist['molecule'][wfilter] == 'T'
            num_molecules = len(np.where(molecules)[0])

            # Append microturbulence, solar abundances and metallicity
            moog_atmosphere = open(atmosphere_filename, "a")
            atom_abundances = abundances[np.logical_and(abundances['code'] > 1, abundances['code'] <= 92)] # Don't update hydrogen or helium abundances
            moog_atmosphere.write("  %.2f\n" % (microturbulence_vel))
            moog_atmosphere.write("NATOMS=   %i %.2f\n" % (len(atom_abundances), MH))
            for atom_abundance in atom_abundances:
                abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
                moog_atmosphere.write("%i  %.2f\n" % (atom_abundance['code'], abund))
            if num_molecules > 0:
                unique_molecules = np.unique(linelist['spectrum_moog_species'][wfilter][molecules])
                moog_atmosphere.write("NMOL      %i\n" % (len(unique_molecules)))
                for specie in unique_molecules:
                    moog_atmosphere.write("  %s\n" % (specie))
                #moog_atmosphere.write("NMOL      22\n")
                #moog_atmosphere.write("  101.0   106.0   107.0   108.0   112.0  126.0\n")
                #moog_atmosphere.write("  606.0   607.0   608.0\n")
                #moog_atmosphere.write("  707.0   708.0\n")
                #moog_atmosphere.write("  808.0   812.0   822.0\n")
                #moog_atmosphere.write("  10108.0 60808.0\n")
                #moog_atmosphere.write("  6.1     7.1     8.1   12.1  22.1  26.1")
            moog_atmosphere.close()

            # Add hydrogen lines
            # Provide some margin or near-by deep lines might be omitted
            margin = 2. # 2 nm
            wfilter = np.logical_and(hydrogen_lines['wave_A'] >= (wave_base-margin)*10., hydrogen_lines['wave_A'] <= (wave_top+margin)*10.)
            selected_hydrogen_lines = hydrogen_lines[wfilter]
            if len(selected_hydrogen_lines) > 40:
                # TODO: Find a work around to this
                raise Exception("Wavelength range too big, it includes too many hydrogen lines")
            out = open(stronglinelist_file, "w")
            for line in selected_hydrogen_lines:
                out.write("%10.3f%10s%10.3f%10.3f%10s%10s%10s%10s\n" \
                        % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_eV'], line['loggf'], "", "", "", ""))
            out.close()

            par_file = open(tmp_execution_dir + "/batch.par", "w")
            par_file.write("synth\n")
            par_file.write("standard_out moog.std\n")
            par_file.write("summary_out  moog.sum\n")
            par_file.write("smoothed_out moog.spec\n")
            par_file.write("model_in     model.in\n")
            par_file.write("lines_in     lines.in\n")
            par_file.write("stronglines_in stronglines.in\n")
            par_file.write("atmosphere  0\n") # controls the output of atmosphere quantities
            #if num_molecules > 0:
                #par_file.write("molecules   1\n") # controls the molecular equilibrium calculations (1 = do molecular equilibrium but do not print results)
            #else:
                #par_file.write("molecules   0\n")
            # molecules should be always on or moog does not run
            par_file.write("molecules   1\n") # controls the molecular equilibrium calculations (1 = do molecular equilibrium but do not print results)
            par_file.write("lines       1\n") # controls the output of line data (print out standard information about the input line list)
            par_file.write("strong      1\n")
            par_file.write("flux/int    0\n") # choses integrated flux or central intensity (0 = integrated flux calculations)
            par_file.write("damping     1\n") # Use Barklem.dat file with waals, alpha and rad damping coeff, if not found then do like damping = 0
            par_file.write("freeform    0\n") # Linelist format of 7 columns with numbers %10.3f and comment %10s
            par_file.write("plot        3\n")
            par_file.write("abundances  0  1\n")
            par_file.write("isotopes    0  1\n")
            par_file.write("synlimits\n")
            wave_range_of_line_influence = 20.0 # Amstrong
            par_file.write(" %.2f %.2f %.2f %.1f\n" % (wave_base*10., wave_top*10., wave_step*10., wave_range_of_line_influence))
            par_file.write("obspectrum    5\n")
            #par_file.write("plotpars    0\n")
            par_file.write("plotpars    1\n")
            par_file.write(" %.2f %.2f 0.0 1.10\n" % (wave_base*10., wave_top*10.)) # Ploting limits
            par_file.write(" 0.0     0.0    0.00  1.0\n") # vshift       lamshift    obsadd    obsmult
            par_file.write(" r       0.0  0.0  0.0  0.00  0.0") #  smooth-type  FWHM-Gauss  vsini     LimbDarkeningCoeff    FWHM-Macro     FWHM-Loren
            par_file.close()

            previous_cwd = os.getcwd()
            os.chdir(tmp_execution_dir)

            command = moog_executable
            command_input = ""

            # If timeout command exists in PATH, then use it to control moog execution time
            if which("timeout") is not None:
                command = "timeout %i " % (timeout) + command

            # MOOG clears the screen, so it is better to always don't do verbose this part
            # also the information printed is not specially useful
            #if verbose == 1:
                #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
            #else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the moog synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            os.chdir(previous_cwd)

            try:
                data = np.loadtxt(tmp_execution_dir+"/moog.spec", skiprows=2)
            except:
                print out
                sys.stdout.flush()
                raise Exception("Synthesis failed!")
            #synth_waveobs_tmp = np.linspace(wave_base, wave_top, len(synth_fluxes_tmp)) # Not exactly identical to turbospectrum wavelengths
            synth_waveobs_tmp = data[:,0] / 10. # Armstrong to nm
            synth_waveobs = np.hstack((synth_waveobs, synth_waveobs_tmp))

            synth_fluxes_tmp = data[:,1]
            if len(synth_fluxes_tmp) > 1 and np.isnan(synth_fluxes_tmp[-1]):
                synth_fluxes_tmp[-1] = synth_fluxes_tmp[-2] # Turbospectrum bug with gfortran, last flux is always NaN
            synth_fluxes = np.hstack((synth_fluxes, synth_fluxes_tmp))

            shutil.rmtree(tmp_execution_dir)

    segments = None
    vrad = (0,)
    synth_fluxes = apply_post_fundamental_effects(synth_waveobs, synth_fluxes, segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)

    synth_spectrum = create_spectrum_structure(synth_waveobs, synth_fluxes)
    synth_spectrum.sort(order=['waveobs'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(synth_spectrum['waveobs'], waveobs):
        synth_spectrum = resample_spectrum(synth_spectrum, waveobs, method="linear", zero_edges=True)

    return synth_spectrum['flux']


def __synthe_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, molecules_files=None, regions=None, tmp_dir=None, timeout=1800):
    return __synthe_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, tmp_dir=tmp_dir, timeout=timeout)

def __synthe_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, molecules_files=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, tmp_dir=None, timeout=1800):
    if not is_synthe_support_enabled():
        raise Exception("Synthe support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    atmos_dir = ispec_dir + "/synthesizer/atmos/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        xnfpelsyn_executable = atmos_dir + "bin.amd64/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.amd64/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.amd64/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.amd64/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.amd64/synthe.exe"
        spectrv_executable = atmos_dir + "bin.amd64/spectrv.exe"
        rotate_executable = atmos_dir + "bin.amd64/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.amd64/syntoascanga.exe"
    else:
        logging.warn("*************************************************")
        logging.warn("Synthe does not work properly in 32 bits systems!")
        logging.warn("*************************************************")
        xnfpelsyn_executable = atmos_dir + "bin.ia32/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.ia32/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.ia32/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.ia32/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.ia32/synthe.exe"
        spectrv_executable = atmos_dir + "bin.ia32/spectrv.exe"
        rotate_executable = atmos_dir + "bin.ia32/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.ia32/syntoascanga.exe"
    atmos_molecules = atmos_dir + "lines/molecules.dat"
    atmos_helium = atmos_dir + "lines/he1tables.dat"
    atmos_continua = atmos_dir + "lines/continua.dat"


    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top
    else:
        global_wave_base = np.min(regions['wave_base'])
        global_wave_top = np.max(regions['wave_top'])

    # TODO: decide how to stablish wave_step
    waveobs = waveobs.copy()
    waveobs.sort()
    #wave_step = waveobs[1] - waveobs[0] # 0.001
    wave_step = np.max((0.001, np.min(waveobs[1:] - waveobs[:-1])))

    # Limit linelist
    linelist = __filter_linelist(linelist, regions)


    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        abundances = abundances.copy()
        for fixed_abundance in fixed_abundances:
            index = np.where(abundances['code'] == fixed_abundance['code'])[0]
            abundances['Abund'][index] = fixed_abundance['Abund'] - MH



    synth_fluxes = []
    synth_waveobs = []
    for i, region in enumerate(regions):
        # It is better not to synthesize in a single run a big chunk of wavelength so
        # we split the computation in several pieces
        max_segment = 100. # nm
        if (region['wave_top'] - region['wave_base'])/wave_step > max_segment/wave_step:
            segment_wave_base = np.arange(region['wave_base'], region['wave_top'], max_segment)
            segments = np.recarray((len(segment_wave_base),),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'] = segment_wave_base
            segments['wave_top'] = segment_wave_base + max_segment - wave_step
            segments['wave_top'][-1] = region['wave_top'] # Last segment should not over pass the original global limits
        else:
            segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'][0] = region['wave_base']
            segments['wave_top'][0] = region['wave_top']

        for segment in segments:
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']

            tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
            os.symlink(atmos_molecules, tmp_execution_dir+"/fort.2")
            os.symlink(atmos_helium, tmp_execution_dir+"/fort.18")
            os.symlink(atmos_continua, tmp_execution_dir+"/fort.17")


            previous_cwd = os.getcwd()
            os.chdir(tmp_execution_dir)

            # XNFPELSYN pretabulates continuum opacities and number densities for
            # different chemical elements and writes them to fort.10
            command = xnfpelsyn_executable

            command_input = "SURFACE INTENSI 17 1.,.9,.8,.7,.6,.5,.4,.3,.25,.2,.15,.125,.1,.075,.05,.025,.01\n"
            #command_input = "SURFACE FLUX\n"
            command_input += "ITERATIONS 1 PRINT 2 PUNCH 2\n"
            command_input += "CORRECTION OFF\n"
            command_input += "PRESSURE OFF\n"
            command_input += "READ MOLECULES\n"
            command_input += "MOLECULES ON\n"
            command_input += "TEFF   %.0f  GRAVITY %.5f LTE\n" % (teff, logg)
            command_input += "TITLE  ISPEC\n"
            command_input += " OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0\n"
            mixing_length_param = 1.25
            command_input += " CONVECTION ON   %.2f TURBULENCE OFF  0.00  0.00  0.00  0.00\n" % (mixing_length_param)
            abundance_scale = 10**MH
            ## Fraction in number of the total atoms from WIDTH example:
            #   hydrogen_number_atom_fraction = 0.92080
            #   helium_number_atom_fraction = 0.07837
            # Fraction in mass from MARCS model (X=Hydrogen, Y=Helium, Z=Metals):
            #   0.74732 0.25260 7.81E-05 are X, Y and Z, 12C/13C=89 (=solar)
            # Transform to number fraction:
            #   Y = 0.25260 / (4-3*0.25260) = 0.07791
            #   X = 1 - Y = 0.92209
            #hydrogen_number_atom_fraction = 0.92209
            #helium_number_atom_fraction = 0.07791
            hydrogen_number_atom_fraction = 0.92080 # It does not seem to have any effect on synthesis
            helium_number_atom_fraction = 0.07837   # It does not seem to have any effect on synthesis

            command_input += "ABUNDANCE SCALE   %.5f ABUNDANCE CHANGE 1 %.5f 2 %.5f\n" % (abundance_scale, hydrogen_number_atom_fraction, helium_number_atom_fraction)
            # command_input += " ABUNDANCE CHANGE  3 -10.99  4 -10.66  5  -9.34  6  -3.65  7  -4.26  8  -3.38\n"
            # command_input += " ABUNDANCE CHANGE  9  -7.48 10  -4.20 11  -5.87 12  -4.51 13  -5.67 14  -4.53\n"
            atom_abundances = abundances[np.logical_and(abundances['code'] > 2, abundances['code'] <= 92)]
            num_added_abundances = 0
            for atom_abundance in atom_abundances:
                # abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
                #command_input += " ABUNDANCE CHANGE  %i  %.2f\n" % (atom_abundance['code'], abund)
                if num_added_abundances == 0:
                    command_input += " ABUNDANCE CHANGE"

                abund = atom_abundance['Abund']
                command_input += " %2i %6.2f" % (atom_abundance['code'], abund)
                num_added_abundances += 1

                if num_added_abundances == 6:
                    command_input += "\n"
                    num_added_abundances = 0
            command_input += " ABUNDANCE CHANGE 93 -20.00 94 -20.00 95 -20.00 96 -20.00 97 -20.00 98 -20.00    \n"
            command_input += " ABUNDANCE CHANGE 99 -20.00                                                      \n"


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
            model = command_input

            # If timeout command exists in PATH, then use it to control synthe execution time
            which_timeout = which("timeout")
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            # Never verbose because WIDTH's output is printed on stdout (not saved on a file)
            # and it is needed by iSpec
            #if verbose == 1:
                #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
            #else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode
            #print out

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            # synbeg: Reads the fundamental parameters of the process: wavelength range,
            # resolution of the final spectrum before degrading, whether unclassified
            # (pre-dicted) lines should be included, etc., and writes them to fort.93
            command = synbeg_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            # microturbulence is added by summing the squares to the microturbulence indicated in the mode atmosphere, which we have forced to zero before
            turbv = microturbulence_vel
            # NOTE: For the wavelength range, it does not work as other codes... it will ignore any lines outside the provided range
            # even if they are strong near-by lines that will affect the region of interest. So we increase the region to synthesize
            # and we will cut later.
            # Provide some margin or near-by deep lines might be omitted
            margin = 2. # 2 nm
            #command_input =  "AIR        %-9.1f %-9.1f 600000. %9.2f    0     30    .0001     1    0\n" % (wave_base-margin, wave_top+margin, turbv)
            command_input =  "AIR        %-9.1f %-9.1f 600000. %9.2f    0     10    .001      0    0\n" % (wave_base-margin, wave_top+margin, turbv)
            command_input += "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD\n"
            #AIR indicates that the wavelengths are in AIR. VAC would provide vacuum wavelengths
            #WLBEG and WLEND are the starting and ending points of the synthesis, in nanometers
            #RESOLU is the resolution at which the calculation is performed. Practically, SYNTHE calculates the transfer through the atmosphere at wavelength intervals with such spacing. Of course, reducing the resolution will lead to a faster calculation, but also to a poorer sampling of the radiative transfer through the atmosphere. We thus suggest not to go below a resolution of 100000. This value is adequate for comparison with high resolution observed spectra.
            #TURBV is the microturbulence we want SYNTHE to add to the one in the atmosphere model. Since microturbulence is added by summing the squares, and we have a VTURB=1 model, we need to add 1.67 to obtain the final 1.95 km/s.
            #IFNLTE is set to 0 because we want a LTE calculation
            #LINOUT if it is negative, line data are not saved and it speeds up the process
            #CUTOFF is used to keep the weakest transitions out of the output files. With this setting, any absorption subtracting at its center less than 1/10000 of the intensity will be cut off.

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode
            #print out

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")


            if linelist_file is None:
                # Provide some margin or near-by deep lines might be omitted
                margin = 2. # 2 nm
                wfilter = np.logical_and(linelist['wave_nm'] >= wave_base-margin, linelist['wave_nm'] <= wave_top+margin)
                linelist_filename, molecules_filenames = write_atomic_linelist(linelist[wfilter], code="synthe", tmp_dir=tmp_execution_dir)
            else:
                linelist_filename = linelist_file
                molecules_filenames = molecules_files

            ## Compute atomic lines
            os.symlink(linelist_filename, tmp_execution_dir+"/fort.11")
            # rgfallinesnew: adds line information and line opacity data by reading
            # the adequate opacity file from fort.11.
            # The line opacity data are sent to fort.12 (and the line identification
            # to fort.14) if the calculation is LTE; it goes to fort.19 (line
            # identification to fort.20) if not.
            command = rgfallinesnew_executable
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            #print out
            #print "-"*80
            errcode = proc.returncode
            os.remove(tmp_execution_dir+"/fort.11")

            # Compute molecules
            if molecules_filenames is not None:
                for molecules_filename in molecules_filenames:
                    os.symlink(molecules_filename, tmp_execution_dir+"/fort.11")
                    # RMOLEC is the homologue of RGFALLTEST for diatomic molecules
                    # It writes line opacity to fort.12 and line identifications
                    # to fort. 14
                    command = rmolescasc_executable
                    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                    # wait for the process to terminate
                    out, err = proc.communicate()
                    errcode = proc.returncode
                    os.remove(tmp_execution_dir+"/fort.11")
                    #print "Filename:", molecules_filename
                    #print out


            # SYNTHE computes line opacity data based on the fundamental parameters of
            # the model (the ones written by SYNBEG to fort.93) and writes them to fort.93
            command = synthe_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            config_data = "0.0       0.        1.        0.        0.        0.        0.        0.\n"
            config_data += "0.\n"
            config_data += "RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW\n"
            config = open(tmp_execution_dir+"/fort.25", "w")
            config.write(config_data)
            config.close()

            # SPECTRV reads from unit 9 the file written by the program SYNTHE and com-
            # putes the continuum opacities and the overall synthetic spectrum
            # (intensities at 17 angles); this spectrum is then written to fort.7
            command = spectrv_executable
            command_input = model

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")



            # ROTATE handles rotational broadening of spectral lines. It also integrates the
            # intensity data given by SPECTRV to compute total flux.
            # The arguments are the number of values of the projected rotational
            # velocity v sin i (first row) and the velocities themselves (second row)
            #
            # Working notes: "It is recommended to always use SURFACE INTENSITY instead of SURFACE FLUX.
            # we can get flux spectra of non-rotating stars with SURFACE INTENSITY and
            # v sin i = 0"
            command = rotate_executable
            intensities_filename = tmp_execution_dir+"/intensities.bin"
            if os.path.exists(tmp_execution_dir+"/fort.7"):
                os.rename(tmp_execution_dir+"/fort.7", intensities_filename)
            else:
                raise Exception("No intensities were calculated, synthesis failed!")
            #os.remove(tmp_execution_dir+"/fort.1")
            os.symlink(intensities_filename, tmp_execution_dir+"/fort.1")
            #command_input = "    1   50\n"
            command_input = "    1\n"
            command_input +="0.\n"

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            flux_filename = tmp_execution_dir+"/flux.bin"
            os.rename(tmp_execution_dir+"/ROT1", flux_filename)
            os.remove(tmp_execution_dir+"/fort.1")
            os.remove(tmp_execution_dir+"/fort.2")

            #os.remove(tmp_execution_dir+"/fort.3")
            os.symlink(flux_filename, tmp_execution_dir+"/fort.1")
            os.symlink(tmp_execution_dir+"/lines.txt", tmp_execution_dir+"/fort.3")
            os.symlink(tmp_execution_dir+"/spectrum.txt", tmp_execution_dir+"/fort.2")
            os.symlink(tmp_execution_dir+"/dump.txt", tmp_execution_dir+"/fort.4")


            command = syntoascanga_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            try:
                data = np.loadtxt(tmp_execution_dir+"/spectrum.txt")
            except:
                #print out
                sys.stdout.flush()
                raise Exception("Synthesis failed!")
            synth_waveobs_tmp = data[:,0] / 10. # Armstrong to nm
            # NOTE: We provided an artificially bigger wavelength range, so that synthe will consider near-by deep lines
            # Now we correct that and we reduce the wavelength range to the correct one:
            wfilter = np.logical_and(synth_waveobs_tmp >= wave_base, synth_waveobs_tmp <= wave_top)
            synth_waveobs = np.hstack((synth_waveobs, synth_waveobs_tmp[wfilter]))

            synth_fluxes_tmp = data[:,3]
            synth_fluxes = np.hstack((synth_fluxes, synth_fluxes_tmp[wfilter]))

            os.chdir(previous_cwd)
            shutil.rmtree(tmp_execution_dir)

    segments = None
    vrad = (0,)
    synth_fluxes = apply_post_fundamental_effects(synth_waveobs, synth_fluxes, segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)

    synth_spectrum = create_spectrum_structure(synth_waveobs, synth_fluxes)
    synth_spectrum.sort(order=['waveobs'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(synth_spectrum['waveobs'], waveobs):
        synth_spectrum = resample_spectrum(synth_spectrum, waveobs, method="linear", zero_edges=True)

    return synth_spectrum['flux']




def __sme_generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, regions=None, timeout=1800):
    return __sme_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose, regions=regions, timeout=timeout, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0)


def __sme_true_generate_spectrum(process_communication_queue, waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None):
    if not is_sme_support_enabled():
        raise Exception("SME support is not enabled")


    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    sme_dir = ispec_dir + "/synthesizer/sme/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        sme = ctypes.CDLL(sme_dir + "/sme_synth.so.linux.x86_64.64")
    else:
        sme = ctypes.CDLL(sme_dir + "/sme_synth.so.linux.x86.32")

    #logging.warn("SME does not support isotope modifications")

    waveobs = waveobs.copy()
    waveobs.sort()
    #wave_step = waveobs[1] - waveobs[0] # 0.001
    wave_step = np.max((0.001, np.min(waveobs[1:] - waveobs[:-1])))

    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top
    else:
        global_wave_base = np.min(regions['wave_base'])
        global_wave_top = np.max(regions['wave_top'])

    # Limit linelist
    linelist = __filter_linelist(linelist, regions)

    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    atom_abundances = abundances[abundances['code'] <= 92]
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        atom_abundances = atom_abundances.copy()
        for fixed_abundance in fixed_abundances:
            index = np.where(atom_abundances['code'] == fixed_abundance['code'])[0]
            # WARNING: Fix the abundance BUT also substract MH because SME will scale it later on
            atom_abundances['Abund'][index] = fixed_abundance['Abund'] - MH

    radius = atmosphere_layers[0][-1]
    nvalues = len(atmosphere_layers[0])
    if nvalues == 11 and radius > 2.0: # Compare to 2.0 instead of 1.0 to avoid floating point imprecisions
        logging.info("Spherical model atmosphere with radius %.2e cm" % (radius))
        spherical_model = True
    else:
        spherical_model = False

    #---------------------------------------------------------------------------
    # *.- Entry.SMELibraryVersion
    #msg = _sme_librrayversion(sme)

    #---------------------------------------------------------------------------
    # 0.- Entry.InputLineList (passlinelist.pro)
    if verbose == 1:
        logging.info("SME InputLineList")
    msg = _sme_inputlinelist(sme, linelist)
    if msg != "''":
        logging.warn(msg)

    #---------------------------------------------------------------------------
    # 1.- Entry.InputModel (passmodel.pro)
    if verbose == 1:
        logging.info("SME InputModel")
    msg = _sme_inputmodel(sme, teff, logg, MH, microturbulence_vel, atmosphere_layers, atom_abundances, spherical_model)
    if msg != "''":
        logging.warn(msg)

    #---------------------------------------------------------------------------
    # 2.- Entry.InputNLTE

    #---------------------------------------------------------------------------
    # 3.- Entry.InputAbund (passabund.pro)
    if verbose == 1:
        logging.info("SME InputAbund")
    msg = _sme_inputabund(sme, atom_abundances, MH)
    if msg != "''":
        logging.warn(msg)

    #---------------------------------------------------------------------------
    # 4.- Entry.Ionization
    if verbose == 1:
        logging.info("SME Ionization")
    msg = _sme_ionization(sme)
    if msg != "''":
        logging.warn(msg)

    #---------------------------------------------------------------------------
    # 5.- Entry.SetVWscale
    if verbose == 1:
        logging.info("SME SetVWscale")
    msg = _sme_setvwscale(sme)
    if msg != "''":
        logging.warn(msg)


    synth_fluxes = []
    synth_waveobs = []
    for i, region in enumerate(regions):
        # It is better not to synthesize in a single run a big chunk of wavelength so
        # we split the computation in several pieces
        max_segment = 100. # nm
        if (region['wave_top'] - region['wave_base'])/wave_step > max_segment/wave_step:
            segment_wave_base = np.arange(region['wave_base'], region['wave_top'], max_segment)
            segments = np.recarray((len(segment_wave_base),),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'] = segment_wave_base
            segments['wave_top'] = segment_wave_base + max_segment - wave_step
            segments['wave_top'][-1] = region['wave_top'] # Last segment should not over pass the original global limits
        else:
            segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'][0] = region['wave_base']
            segments['wave_top'][0] = region['wave_top']

        for segment in segments:
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']
            if verbose == 1:
                logging.info("Segment %.2f - %.2f" % (wave_base, wave_top))


            #---------------------------------------------------------------------------
            # 6.- Entry.InputWaveRange (passwaverange.pro)
            if verbose == 1:
                logging.info("SME InputWaveRange")
            msg = _sme_inputwaverange(sme, wave_base*10., wave_top*10.)
            if msg != "''":
                logging.warn(msg)

            #---------------------------------------------------------------------------
            # 7.- Entry.Opacity
            if verbose == 1:
                logging.info("SME Opacity")
            msg = _sme_opacity(sme)
            if msg != "''":
                logging.warn(msg)

            #---------------------------------------------------------------------------
            # 8.- Entry.Transf
            if verbose == 1:
                logging.info("SME Transf")
            first_execution = i == 0
            #nwmax = int((wave_top - wave_base) / wave_step) * 2
            nwmax = 200000
            synth_waveobs_tmp, synth_fluxes_tmp = _sme_transf(sme, sme_dir, nwmax, keep_lineop=not first_execution)
            if msg != "''":
                logging.warn(msg)

            synth_waveobs = np.hstack((synth_waveobs, synth_waveobs_tmp))
            synth_fluxes = np.hstack((synth_fluxes, synth_fluxes_tmp))

    segments = None
    vrad = (0,)
    synth_fluxes = apply_post_fundamental_effects(synth_waveobs, synth_fluxes, segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)

    synth_spectrum = create_spectrum_structure(synth_waveobs/10., synth_fluxes)
    synth_spectrum.sort(order=['waveobs'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(synth_spectrum['waveobs'], waveobs):
        synth_spectrum = resample_spectrum(synth_spectrum, waveobs, method="linear", zero_edges=True)

    process_communication_queue.put(synth_spectrum['flux'])



def __sme_generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, timeout=1800):
    if not is_sme_support_enabled():
        raise Exception("SME support is not enabled")

    fluxes = np.zeros(len(waveobs))

    # It is better to run SME in a separate process since sometimes it
    # aborts the full process because unknown reasons
    process_communication_queue = Queue()

    p = Process(target=__sme_true_generate_spectrum, args=(process_communication_queue, waveobs, atmosphere_layers, teff, logg, MH, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel), kwargs={'regions': regions, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'verbose': verbose})
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
            #elif gui_queue is not None:
                ## GUI update
                ## It allows communications between process in order to update the GUI progress bar
                #gui_queue.put(data)
            process_communication_queue.task_done()
        except Empty:
            # No results, continue waiting
            pass
        num_seconds += 1
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the SME synthesis process.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("SME has failed.")
        p.terminate()
    else:
        p.join()

    return fluxes


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Functions for vsini and vmac broadening
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

def __lsf_rotate(deltav,vsini,epsilon=0.6):
    # Based on lsf_rotate.pro:
    #  http://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro
    #
    # Adapted from rotin3.f in the SYNSPEC software of Hubeny & Lanz
    # http://nova.astro.umd.edu/index.html    Also see Eq. 17.12 in
    # "The Observation and Analysis of Stellar Photospheres" by D. Gray (1992)
    e1 = 2.0*(1.0 - epsilon)
    e2 = np.pi*epsilon/2.0
    e3 = np.pi*(1.0 - epsilon/3.0)

    npts = np.ceil(2*vsini/deltav)
    if npts % 2 == 0:
        npts += 1
    nwid = np.floor(npts/2)
    x = np.arange(npts) - nwid
    x = x*deltav/vsini
    x1 = np.abs(1.0 - x**2)

    velgrid = x*vsini
    return velgrid, (e1*np.sqrt(x1) + e2*x1)/e3


def __vmac_broadening(flux, velocity_step, vmac):
    """
        velocity_step: fluxes should correspond to a spectrum homogeneously sampled in velocity space
                    with a fixed velocity step [km/s]
        vmac   : macroturbulence velocity [km/s]

        Based on SME's rtint
        It uses radial-tangential instead of isotropic Gaussian macroturbulence.
    """
    if vmac is not None and vmac > 0:
        # mu represent angles that divide the star into equal area annuli,
        # ordered from disk center (mu=1) to the limb (mu=0).
        # But since we don't have intensity profiles at various viewing (mu) angles
        # at this point, we just take a middle point:
        m = 0.5
        # Calc projected simga for radial and tangential velocity distributions.
        sigma = vmac/np.sqrt(2.0) / velocity_step
        sigr = sigma * m
        sigt = sigma * np.sqrt(1.0 - m**2.)
        # Figure out how many points to use in macroturbulence kernel
        nmk = max(min(round(sigma*10), (len(flux)-3)/2), 3)
        # Construct radial macroturbulence kernel w/ sigma of mu*vmac/sqrt(2)
        if sigr > 0:
            xarg = (np.arange(2*nmk+1)-nmk) / sigr   # exponential arg
            #mrkern = np.exp(max((-0.5*(xarg**2)),-20.0))
            mrkern = np.exp(-0.5*(xarg**2))
            mrkern = mrkern/mrkern.sum()
        else:
            mrkern = np.zeros(2*nmk+1)
            mrkern[nmk] = 1.0    #delta function

        # Construct tangential kernel w/ sigma of sqrt(1-mu**2)*vmac/sqrt(2.)
        if sigt > 0:
            xarg = (np.arange(2*nmk+1)-nmk) /sigt
            mtkern = np.exp(-0.5*(xarg**2))
            mtkern = mtkern/mtkern.sum()
        else:
            mtkern = np.zeros(2*nmk+1)
            mtkern[nmk] = 1.0

        ## Sum the radial and tangential components, weighted by surface area
        area_r = 0.5
        area_t = 0.5
        mkern = area_r*mrkern + area_t*mtkern

        # Convolve the flux with the kernel
        flux_conv = 1 - fftconvolve(1-flux, mkern, mode='same') # Fastest
        #import scipy
        #flux_conv = scipy.convolve(flux, mkern, mode='same') # Equivalent but slower

        return flux_conv
    else:
        return flux

def __vsini_broadening_limbdarkening(flux, velocity_step, vsini, epsilon):
    """
        velocity_step: fluxes should correspond to a spectrum homogeneously sampled in velocity space
                    with a fixed velocity step [km/s]
        vsini   : rotation velocity [km/s]
        epsilon : numeric scalar giving the limb-darkening coefficient,
               default = 0.6 which is typical for  photospheric lines.

        Based on lsf_rotate.pro:
        http://idlastro.gsfc.nasa.gov/ftp/pro/astro/lsf_rotate.pro

        Adapted from rotin3.f in the SYNSPEC software of Hubeny & Lanz
        http://nova.astro.umd.edu/index.html    Also see Eq. 17.12 in
        "The Observation and Analysis of Stellar Photospheres" by D. Gray (1992)
    """
    if vsini is not None and vsini > 0:
        if epsilon is None:
            epsilon = 0.
        kernel_x, kernel_y = __lsf_rotate(velocity_step, vsini, epsilon=epsilon)
        kernel_y /= kernel_y.sum()

        #-- convolve the flux with the kernel
        flux_conv = 1 - fftconvolve(1-flux, kernel_y, mode='same') # Fastest
        #import scipy
        #flux_conv = 1 - scipy.convolve(1-flux, kernel_y, mode='same') # Equivalent but slower
        return flux_conv
    else:
        return flux

def __determine_velocity_step(spectrum):
    # Determine step size for a new model wavelength scale, which must be uniform
    # in velocity to facilitate convolution with broadening kernels. The uniform
    # step size is the largest of:
    wave_base = spectrum['waveobs'][0]
    wave_top = spectrum['waveobs'][-1]
    wmid = (wave_top + wave_base) / 2. # midpoint
    wspan = wave_top - wave_base # width
    # Light speed in vacuum
    #c = 299792458.0 # m/s
    c = 299792.4580 # km/s


    # [1] smallest wavelength step considering the wavelength sampling
    wave_diff = spectrum['waveobs'][1:] - spectrum['waveobs'][:-1]
    min_wave_step = np.min(wave_diff)
    min_wave_step_index = np.argmin(wave_diff)
    vstep1 = min_wave_step / (spectrum['waveobs'][min_wave_step_index] * c)

    # [2] 10% the mean dispersion
    vstep2 = 0.1 * wspan / len(spectrum) / (wmid * c)

    # [3] 0.05 km/s, which is 1% the width of solar line profiles
    vstep3 = 0.05e0

    # Select the largest between 1, 2 and 3:
    velocity_step = np.max((vstep1, vstep2, vstep3))

    return velocity_step
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
