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
import numpy as np
from atmospheres import *
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
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
    if filename is not None:
        out = open(filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    #out.write("\n".join([" ".join(map(str, (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['ew'], line['element']))) for line in linemasks]))
    for line in linemasks:
        # The format is different depending on the broadening parameters
        if line['transition type'] == "AO":
            # O'Mara
            text = "%.5f %s %i %i %f %.1f %s %.4f %f %s" % (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['ew'], line['element'])
        elif line['transition type'] == "GA":
            # Rad, Stark and Waals
            text = "%.5f %s %i %i %f %.1f %s %.4f %.4f %.4f %f %s" % (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['ew'], line['element'])
        else:
            # For i.e. line['transition type'] == "99"
            # Let SPECTRUM calculate them
            text = "%.5f %s %i %i %f %.1f %s %f %s" % (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'],line['log(gf)'], line['fudge factor'], line['transition type'], line['ew'], line['element'])
        out.write(text + "\n")

    out.close()
    return out.name

def __get_element_specie(element_name, chemical_elements, molecules=None):
    """
    Convert element names type "Fe 1" or "Fe 2" to species code form by the atomic number + "." + ionization state
    (i.e. "Fe 2" -> "26.1").
    Or if the elenemt name does not include the ionization state, the returned code correspond only to the atomic_number
    (i.e. "Fe" -> "26").
    Raises an exception if not found.
    """
    element = element_name.split() # Do not specify " " to avoid problems with elements with double spaces like "V  1"

    # Element not present or with a bad format, skip
    if element_name == "":
        raise Exception("Empty element name!")

    symbol = element[0]
    try:
        element.remove('') # Make sure there are not additional spaces between the symbol and the ionization state
        element.remove('')
        element.remove('')
    except ValueError as e:
        pass


    tfilter = (chemical_elements['symbol'] == symbol)
    if len(chemical_elements[tfilter]["atomic_num"]) == 0:
        if molecules is not None:
            # Symbol not found, maybe it is a molecule
            mfilter = (molecules['symbol'] == symbol)
            if len(molecules[mfilter]["atomic_num"]) == 0:
                raise Exception("Unkown '%s' element!" % symbol)
            else:
                specie = str(molecules[mfilter]["atomic_num"][0])
        else:
            raise Exception("Unkown '%s' element!" % symbol)
    else:
        specie = str(chemical_elements[tfilter]["atomic_num"][0])

    if len(element) == 2:
        ionization = str(int(element[1]) - 1)
        specie = specie + "." + ionization

    return specie


def create_free_abundances_structure(free_abundance_elements, chemical_elements, solar_abundances):
    """
    Create the needed structure to determine elemental abundances for a given
    list of elements (i.e. ["Fe"] or ["Fe", "Mg", "Ca"]).
    """
    free_abundances = np.recarray((len(free_abundance_elements), ), dtype=[('code', int),('Abund', float), ('element', '|S30')])
    for element_name in free_abundance_elements:
        specie = __get_element_specie(element_name, chemical_elements, molecules=None)
        if "." in specie:
            raise Exception("Bad format '%s'" % element_name)
        free_abundances['code'] = int(specie)
        free_abundances['Abund'] = solar_abundances['Abund'][solar_abundances['code'] == int(specie)]
        free_abundances['element'] = element_name
    return free_abundances

def write_fixed_abundances(fixed_abundances, filename=None):
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
    if filename is not None:
        out = open(filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("TOTAL\n")
    out.write("\n".join([" ".join(map(str, (line['code'], line['Abund']))) for line in fixed_abundances]))
    out.close()
    return out.name


def read_solar_abundances(abundances_filename):
    """
    Load a SPECTRUM format abundances file for spectral synthesis

    Abund field should be in terms of number densities relative to the total number density
    (log(N/Ntot) scale).
    """
    abundances = np.array([tuple(line.rstrip('\r\n').split()) for line in open(abundances_filename,)][1:], \
                            dtype=[('code', int), ('Abund', '<f8'), ('Amass', '<f8'), \
                            ('I1/D0', '<f8'), ('I2/rdmass', '<f8'), ('I3', '<f8'), \
                            ('I4', '<f8'), ('maxcharge', int)])
    return abundances

def write_solar_abundances(abundances, abundances_filename=None):
    """
    Saves a SPECTRUM abundances file for spectral synthesis

    Abund field should be in terms of number densities relative to the total number density
    (log(N/Ntot) scale).

    If filename is not specified, a temporary file is created and the name is returned.
    """
    if abundances_filename is not None:
        out = open(abundances_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("code   Abund   Amass   I1/D0   I2/rdmass     I3         I4     maxcharge\n")
    out.write("\n".join(["  ".join(map(str, (line['code'], line['Abund'], line['Amass'], line['I1/D0'], line['I2/rdmass'], line['I3'], line['I4'], line['maxcharge']))) for line in abundances]))
    out.close()
    return out.name


def determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, verbose=0, gui_queue=None, timeout=900):
    """
    Determine abundances from equivalent widths (linemasks previously fitted and
    cross-matched with an atomic linelist).

    Returns:

    - The abundance on the scale native to SPECTRUM
    - The abundance in he normal scale, in which the logarithmic abundance of hydrogen is equal to 12.0.
    - The abundance relative to the unscaled abundances (specified in "abundances").
      If "abundances" contain solar abundances, these values represent
      the quantity [X/H] where X is the species in question.
    """

    linemasks_file = write_abundance_lines(linemasks)
    atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
    abundances_file = write_solar_abundances(abundances)
    num_measures = len(linemasks)
    nlayers = len(atmosphere_layers)

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__determine_abundances, args=(process_communication_queue, atmosphere_layers_file, linemasks_file, num_measures, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'verbose': verbose})
    p.start()
    # Default values
    spec_abund = np.zeros(num_measures)
    normal_abund = np.zeros(num_measures)
    x_over_h = np.zeros(num_measures)
    x_over_fe = np.zeros(num_measures)
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) == tuple:
                # Results received!
                spec_abund, normal_abund, x_over_h = data
                x_over_fe = x_over_h - MH
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

    os.remove(atmosphere_layers_file)
    os.remove(linemasks_file)
    os.remove(abundances_file)


    return spec_abund, normal_abund, x_over_h, x_over_fe

def __enqueue_progress(process_communication_queue, v):
    process_communication_queue.put(("self.update_progress(%i)" % v))
    process_communication_queue.join()

def __determine_abundances(process_communication_queue, atmosphere_model_file, linelist_file, num_measures, abundances_file, microturbulence_vel = 2.0, nlayers=56, verbose=0):
    """
    Determine abundances from equivalent widths (linemasks previously fitted and
    cross-matched with an atomic linelist).

    Returns:

    - The abundance on the scale native to SPECTRUM
    - The abundance in he normal scale, in which the logarithmic abundance of hydrogen is equal to 12.0.
    - The abundance relative to the unscaled abundances (specified in "abundances").
      If "abundances" contain solar abundances, these values represent
      the quantity [X/H] where X is the species in question.
    """
    import synthesizer
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    abundances = synthesizer.abundances(atmosphere_model_file, linelist_file, num_measures, abundances_file, microturbulence_vel, nlayers, verbose, update_progress_func)
    process_communication_queue.put(abundances)

