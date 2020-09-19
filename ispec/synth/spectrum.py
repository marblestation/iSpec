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
import numpy as np
from multiprocessing import Process
#from multiprocessing import Queue
from multiprocessing import JoinableQueue
from queue import Empty
import logging

from ispec.abundances import write_solar_abundances, write_fixed_abundances, enhance_solar_abundances
from ispec.atmospheres import write_atmosphere
from ispec.lines import write_atomic_linelist, write_isotope_data
from ispec.common import is_spectrum_support_enabled
from .effects import _filter_linelist, apply_post_fundamental_effects


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=1800, tmp_dir=None, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None):
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
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=0.0, vsini=0.0, limb_darkening_coeff=0.00, R=0, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions, tmp_dir=tmp_dir)


def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=0., vsini=0., limb_darkening_coeff=0., R=0, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, isotope_file=None, regions=None, tmp_dir=None):
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
    if len(linelist) > 1000000:
        raise Exception("Linelist too big for SPECTRUM: %i (limit 1000000)" % (len(linelist)))
    if fixed_abundances is None:
        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float), ('element', '|U30')])

    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top

    waveobs_mask = _create_waveobs_mask(waveobs, regions)
    # Limit linelist
    linelist = _filter_linelist(linelist, regions)

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


def __spectrum_true_generate_spectrum(process_communication_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, isotope_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence=0., vsini=0., limb_darkening_coeff=0., R=0, nlayers=56, verbose=0):
    """
    Generate synthetic spectrum and apply macroturbulence, rotation (visini), limb darkening coeff and resolution except
    if all those parameters are set to zero, in that case the fundamental synthetic spectrum is returned.
    """
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    import ispec.synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    ## The convolution (R), rotation broadening (vsini) and macroturbulence broadening (vmac),
    ## do not seem to work as expected in the SPECTRUM code, thus we set them to zero and
    ## we use a python implementation
    #fluxes = ispec.synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, nlayers, verbose, update_progress_func)
    if verbose:
        verbose = 1
    else:
        verbose = 0
    fluxes = ispec.synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file.encode('utf-8'), linelist_file.encode('utf-8'), isotope_file.encode('utf-8'), abundances_file.encode('utf-8'), fixed_abundances_file.encode('utf-8'), microturbulence_vel, 0, 0, 0, 0, nlayers, verbose, update_progress_func)

    # Zero values, when convolved, remain zero so we give a very tiny flux to avoid this problem
    fluxes[fluxes <= 0] = 10e-9

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

    import ispec.synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    output_wave, output_code, output_ew, output_depth = ispec.synthesizer.calculate_ew_and_depth(atmosphere_model_file.encode('utf-8'), linelist_file.encode('utf-8'), isotope_file.encode('utf-8'), abundances_file.encode('utf-8'), num_lines, microturbulence_vel, nlayers, start, end, verbose, update_progress_func)
    process_communication_queue.put((output_wave, output_code, output_ew, output_depth))

def __enqueue_progress(process_communication_queue, v):
    process_communication_queue.put(("self.update_progress(%i)" % v))
    process_communication_queue.join()


def calculate_theoretical_ew_and_depth(atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, microturbulence_vel = 2.0, atmosphere_layers_file=None, abundances_file=None, linelist_file=None, isotope_file=None, verbose=0, gui_queue=None, timeout=1800, tmp_dir=None):
    """
    """
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    supported = linelist['spectrum_support'] == "T"

    abundances = enhance_solar_abundances(abundances, alpha)

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



