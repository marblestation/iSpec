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

def write_abundance_lines(linemasks, filename):
    """
    Write line regions file with the following format:
    ::

        4700.158 26.0 34843 56120 -1.031 1.0 GA 8.18 -4.4  -7.35  61.13 Fe 1
        4701.044 26.0 29730 51002 -1.857 1.0 GA 8.36 -5.07 -7.245 53.63 Fe 1
        4704.948 26.0 29730 50984 -1.57  1.0 GA 8.14 -5.26 -7.246 64.01 Fe 1

    By order: wavelength (it should be in Angstrom), species code, lower and upper excitation estate, loggf,
    fudge factor, transition type, broadening parameters (rad, stark, waals),
    equivalent width (it should be in mAngstrom, 1 A = 1000 mA) and element name

    """
    out = open(filename, "w")
    out.write("\n".join([" ".join(map(str, (line['VALD_wave_peak'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['ew'], line['element']))) for line in linemasks]))
    out.close()

def determine_abundances(atmosphere_model_file, linelist_file, num_measures, abundances_file, microturbulence_vel = 2.0, nlayers=56, verbose=0, update_progress_func=None, timeout=600):
    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    result_queue = Queue()

    # TODO: Allow communications between process in order to update the GUI progress bar
    update_progress_func = None

    p = Process(target=__determine_abundances, args=(result_queue, atmosphere_model_file, linelist_file, num_measures, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'verbose': verbose, 'update_progress_func':update_progress_func})
    p.start()
    try:
        abundances = result_queue.get(timeout=timeout)
    except Empty:
        logging.error("Timeout in the abundance determination!")
        abundances = (np.zeros(num_measures), np.zeros(num_measures), np.zeros(num_measures))
        p.terminate()
    else:
        p.join()
    p.join()
    return abundances

def __determine_abundances(result_queue, atmosphere_model_file, linelist_file, num_measures, abundances_file, microturbulence_vel = 2.0, nlayers=56, verbose=0, update_progress_func=None):
    """
    Determine abundances from equivalent widths

    Returns:

    - The abundance on the scale native to SPECTRUM
    - The abundance in he normal scale, in which the logarithmic abundance of hydrogen is equal to 12.0.
    - The abundance relative to the unscaled abundances (content of "abundances_file").
    If "abundances_file" contain solar abundances, this values represent
    the quantity [X/H] where X is the species in question.
    """
    import synthesizer
    abundances = synthesizer.abundances(atmosphere_model_file, linelist_file, num_measures, abundances_file, microturbulence_vel, nlayers, verbose, update_progress_func)
    result_queue.put(abundances)

