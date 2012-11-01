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
from atmospheres import *
from multiprocessing import Process
from multiprocessing import Queue
from Queue import Empty

import log
import logging

def generate_spectrum(waveobs, atmosphere_model_file, linelist_file, abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, nlayers=56, verbose=0, update_progress_func=None, timeout=600):
    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    result_queue = Queue()

    # TODO: Allow communications between process in order to update the GUI progress bar
    update_progress_func = None

    p = Process(target=__generate_spectrum, args=(result_queue, waveobs, atmosphere_model_file, linelist_file, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose, 'update_progress_func':update_progress_func})
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
    return fluxes

def __generate_spectrum(result_queue, waveobs, atmosphere_model_file, linelist_file, abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, nlayers=56, verbose=0, update_progress_func=None):
    """
    Generate synthetic spectrum.
    """
    import synthesizer
    fluxes = synthesizer.spectrum(waveobs, atmosphere_model_file, linelist_file, abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, nlayers, verbose, update_progress_func)
    result_queue.put(fluxes)

