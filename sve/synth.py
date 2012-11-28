#
#    This file is part of Spectra Visual Editor (SVE).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    SVE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SVE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with SVE. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from abundances import *
from atmospheres import *
from multiprocessing import Process
from multiprocessing import Queue
from Queue import Empty

import log
import logging

def read_SPECTRUM_linelist(linelist_filename):
    """
    Load a SPECTRUM linelist for spectral synthesis
    """
    linelist = np.array([tuple(line.rstrip('\r\n').split()) for line in open(linelist_filename,)], \
                            dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), \
                            ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'), \
                            ('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
                            ('note', '|S100')])
    return linelist

def write_SPECTRUM_linelist(linelist, linelist_filename=None):
    """
    Saves a SPECTRUM linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    if linelist_filename != None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("\n".join(["  ".join(map(str, (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'], line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['note']))) for line in linelist]))
    out.close()
    return out.name



def generate_spectrum(waveobs, waveobs_mask, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, update_progress_func=None, timeout=600, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs only
    if waveobs_mask contains the value 1.0 at the same position.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be empty.
    """

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_fixed_abund_file = False
    remove_tmp_linelist_file = False
    if atmosphere_layers_file == None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
        remove_tmp_atm_file = True
    if abundances_file == None:
        abundances_file = write_SPECTRUM_abundances(abundances)
        remove_tmp_abund_file = True
    if fixed_abundances_file == None:
        fixed_abundances_file = write_SPECTRUM_fixed_abundances(fixed_abundances)
        remove_tmp_fixed_abund_file = True
    if linelist_file == None:
        linelist_file = write_SPECTRUM_linelist(linelist)
        remove_tmp_linelist_file = True
    nlayers = len(atmosphere_layers)

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    result_queue = Queue()

    # TODO: Allow communications between process in order to update the GUI progress bar
    update_progress_func = None

    p = Process(target=__generate_spectrum, args=(result_queue, waveobs, waveobs_mask, atmosphere_layers_file, linelist_file, abundances_file, fixed_abundances_file), kwargs={'microturbulence_vel': microturbulence_vel, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose, 'update_progress_func':update_progress_func})
    p.start()
    fluxes = np.zeros(len(waveobs))
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            fluxes = result_queue.get(timeout=1)
        except Empty:
            # No results, continue waiting
            num_seconds += 1
        else:
            # Results received!
            break
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
    return fluxes

def __generate_spectrum(result_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, nlayers=56, verbose=0, update_progress_func=None):
    """
    Generate synthetic spectrum.
    """
    import synthesizer
    fluxes = synthesizer.spectrum(waveobs, waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, nlayers, verbose, update_progress_func)
    result_queue.put(fluxes)

