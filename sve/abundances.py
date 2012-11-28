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

def write_abundance_lines(linemasks, filename=None):
    """
    Write line regions file with the following format:
    ::

        4700.158 26.0 34843 56120 -1.031 1.0 GA 8.18 -4.4  -7.35  61.13 Fe 1
        4701.044 26.0 29730 51002 -1.857 1.0 GA 8.36 -5.07 -7.245 53.63 Fe 1
        4704.948 26.0 29730 50984 -1.57  1.0 GA 8.14 -5.26 -7.246 64.01 Fe 1

    By order: wavelength (it should be in Angstrom), species code, lower and upper excitation estate, loggf,
    fudge factor, transition type, broadening parameters (rad, stark, waals),
    equivalent width (it should be in mAngstrom, 1 A = 1000 mA) and element name

    If filename is not specified, a temporary file is created and the name is returned.

    """
    if filename != None:
        out = open(filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("\n".join([" ".join(map(str, (line['VALD_wave_peak'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['ew'], line['element']))) for line in linemasks]))
    out.close()
    return out.name

def write_SPECTRUM_fixed_abundances(fixed_abundances, filename=None):
    """
    Write a fixed abundances file the following format:
    ::

        TOTAL
        20 -6.38

    Where the first column is the specie code and the second the fixed abundance.
    in terms of number densities relative to the total number density
    (log(N/Ntot) scale)

    If filename is not specified, a temporary file is created and the name is returned.

    """
    if filename != None:
        out = open(filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("TOTAL\n")
    out.write("\n".join([" ".join(map(str, (line['code'], line['Abund']))) for line in fixed_abundances]))
    out.close()
    return out.name


def read_SPECTRUM_abundances(abundances_filename):
    """
    Load a SPECTRUM abundances file for spectral synthesis

    Abund field should be in terms of number densities relative to the total number density
    (log(N/Ntot) scale).
    """
    abundances = np.array([tuple(line.rstrip('\r\n').split()) for line in open(abundances_filename,)][1:], \
                            dtype=[('code', int), ('Abund', '<f8'), ('Amass', '<f8'), \
                            ('I1/D0', '<f8'), ('I2/rdmass', '<f8'), ('I3', '<f8'), \
                            ('I4', '<f8'), ('maxcharge', int)])
    return abundances

def write_SPECTRUM_abundances(abundances, abundances_filename=None):
    """
    Saves a SPECTRUM abundances file for spectral synthesis

    Abund field should be in terms of number densities relative to the total number density
    (log(N/Ntot) scale).

    If filename is not specified, a temporary file is created and the name is returned.
    """
    if abundances_filename != None:
        out = open(abundances_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("code   Abund   Amass   I1/D0   I2/rdmass     I3         I4     maxcharge\n")
    out.write("\n".join(["  ".join(map(str, (line['code'], line['Abund'], line['Amass'], line['I1/D0'], line['I2/rdmass'], line['I3'], line['I4'], line['maxcharge']))) for line in abundances]))
    out.close()
    return out.name


def determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, verbose=0, update_progress_func=None, timeout=600):

    linemasks_file = write_abundance_lines(linemasks)
    atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
    abundances_file = write_SPECTRUM_abundances(abundances)
    num_measures = len(linemasks)
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

    p = Process(target=__determine_abundances, args=(result_queue, atmosphere_layers_file, linemasks_file, num_measures, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'verbose': verbose, 'update_progress_func':update_progress_func})
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

    os.remove(atmosphere_layers_file)
    os.remove(linemasks_file)
    os.remove(abundances_file)

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

