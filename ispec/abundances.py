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
from lines import write_atomic_linelist, van_der_Waals_ABO_to_single_gamma_format
from common import is_turbospectrum_support_enabled, is_spectrum_support_enabled, is_moog_support_enabled, is_width_support_enabled
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
from Queue import Empty
import subprocess
import shutil

import log
import logging
import re



def __spectrum_write_abundance_lines(linemasks, filename=None, tmp_dir=None):
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
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
    with_ew = 'ew' in linemasks.dtype.names
    line_ew = 0
    for line in linemasks:
        if with_ew:
            line_ew = line['ew']
        # The format is different depending on the broadening parameters
        if line['spectrum_transition_type'] == "AO":
            # O'Mara
            text = "%.5f %s %i %i %f %.1f %s %.4f %f %s" % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_cm1'], line['upper_state_cm1'],line['loggf'], line['spectrum_fudge_factor'], line['spectrum_transition_type'], line['rad'], line_ew, line['element'])
        elif line['spectrum_transition_type'] == "GA":
            # Rad, Stark and Waals
            text = "%.5f %s %i %i %f %.1f %s %.4f %.4f %.4f %f %s" % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_cm1'], line['upper_state_cm1'],line['loggf'], line['spectrum_fudge_factor'], line['spectrum_transition_type'], line['rad'], line['stark'], line['waals'], line_ew, line['element'])
        else:
            # For i.e. line['spectrum_transition_type'] == "99"
            # Let SPECTRUM calculate them
            text = "%.5f %s %i %i %f %.1f %s %f %s" % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_cm1'], line['upper_state_cm1'],line['loggf'], line['spectrum_fudge_factor'], line['spectrum_transition_type'], line_ew, line['element'])
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
    for i, element_name in enumerate(free_abundance_elements):
        specie = __get_element_specie(element_name, chemical_elements, molecules=None)
        if "." in specie:
            raise Exception("Bad format '%s'" % element_name)
        free_abundances['code'][i] = int(specie)
        free_abundances['Abund'][i] = solar_abundances['Abund'][solar_abundances['code'] == int(specie)]
        free_abundances['element'][i] = element_name
    return free_abundances

def write_fixed_abundances(fixed_abundances, filename=None, tmp_dir=None):
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
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
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


def determine_abundance_enchancements(MH, scale=None):
    """
    Calculates alpha_enhancement
    from the metallicity (MH). By default the standard composition of MARCS
    model atmospheres is used, which is the most logical for most analysis.
    """
    if scale is None:
        # MARCS standard scale
        # - http://marcs.astro.uu.se/
        scale = [(+1.00,   +0.00),
                 (+0.75,   +0.00),
                 (+0.50,   +0.00),
                 (+0.25,   +0.00),
                 (+0.00,   +0.00),
                 (-0.25,   +0.10),
                 (-0.50,   +0.20),
                 (-0.75,   +0.30),
                 (-1.00,   +0.40),
                 (-1.50,   +0.40),
                 (-2.00,   +0.40),
                 (-2.50,   +0.40),
                 (-3.00,   +0.40),
                 (-4.00,   +0.40),
                 (-5.00,   +0.40)]
        scale = np.array(scale, dtype=[('[Fe/H]', float), ('[a/Fe]', float)])
    scale.sort(order='[Fe/H]') # It should be sorted or interpolation does not work properly
    alpha_enhancement = np.interp(MH, scale['[Fe/H]'], scale['[a/Fe]'])
    return alpha_enhancement



# TODO
#abundances = enhance_solar_abundances(abundances, alpha, MH_compensantion=MH)
def enhance_solar_abundances(abundances, alpha_enhancement, MH_compensantion=0.):
    """
    Scales alpha elements and CNO abundances.
    """
    abundances = abundances.copy()

    #  6|C|Carbon|14|2|6|4|4
    c = abundances['code'] == 6
    #  7|N|Nitrogen|15|2|7|6|8
    n = abundances['code'] == 7
    #  8|O|Oxygen|16|2|8|3|3
    o = abundances['code'] == 8


    # APOGEE alpha element: O, Mg, Si, S, Ca, and Ti
    ## http://www.sdss.org/dr12/irspec/aspcap/
    ## http://www.iac.es/proyecto/ATLAS-APOGEE/
    # Kurucz alpha enchanced: O, Ne, Mg, Si, S, Ar, Ca, and Ti
    ## http://kurucz.harvard.edu/grids.html

    # 10|Ne|Neon|18|2|10|5|5
    # 12|Mg|Magnesium|2|3|12|7|6
    alpha = np.logical_or(abundances['code'] == 10, abundances['code'] == 12)
    # 14|Si|Silicon|14|3|14|8|7
    alpha = np.logical_or(alpha, abundances['code'] == 14)
    # 16|S|Sulfur|16|3|16|10|9
    alpha = np.logical_or(alpha, abundances['code'] == 16)
    # 18|Ar|Argon|18|3|18|14|11
    alpha = np.logical_or(alpha, abundances['code'] == 18)
    # 20|Ca|Calcium|2|4|20|12|12
    alpha = np.logical_or(alpha, abundances['code'] == 20)
    # 22|Ti|Titanium|4|4|22|22|19
    alpha = np.logical_or(alpha, abundances['code'] == 22)

    ## MH_compensantion is specified when it is known that abundances will
    ## be scaled by MH after the alpha enhancement (thus it is necessary
    ## to cancel that out)
    #abundances['Abund'][c] += c_enhancement - MH_compensantion
    #abundances['Abund'][n] += n_enhancement - MH_compensantion
    #abundances['Abund'][o] += o_enhancement - MH_compensantion
    abundances['Abund'][o] += alpha_enhancement - MH_compensantion # MARCS standard enhances Oxygen as well as alpha elements at the same rate
    abundances['Abund'][alpha] += alpha_enhancement - MH_compensantion

    return abundances


def write_solar_abundances(abundances, abundances_filename=None, tmp_dir=None):
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
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
    out.write("code   Abund   Amass   I1/D0   I2/rdmass     I3         I4     maxcharge\n")
    out.write("\n".join(["  ".join(map(str, (line['code'], line['Abund'], line['Amass'], line['I1/D0'], line['I2/rdmass'], line['I3'], line['I4'], line['maxcharge']))) for line in abundances]))
    out.close()
    return out.name


def determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, gui_queue=None, timeout=1800, isotopes=None, code="spectrum", tmp_dir=None):
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'width']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    abundances = enhance_solar_abundances(abundances, alpha)

    if code == "turbospectrum":
        return __turbospectrum_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, tmp_dir=tmp_dir)
    elif code == "moog":
        #return __moog_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, tmp_dir=tmp_dir)
        success = False
        bad = linemasks['wave_A'] < 0 # All to false
        # MOOG can fail for some lines and it stops the execution, automatically detect that
        # filter the line and re-execute
        while not success and np.any(~bad):
            try:
                if ignore is not None:
                    filtered_ignore = ignore[~bad]
                else:
                    filtered_ignore = None
                spec_abund, absolute_abund, x_over_h, x_over_fe = __moog_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks[~bad], isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=filtered_ignore, verbose=verbose, tmp_dir=tmp_dir)
            except Exception, e:
                # MOOG ERROR: CANNOT DECIDE ON LINE WAVELENGTH STEP SIZE FOR   5158.62   I QUIT!
                if "CANNOT DECIDE ON LINE WAVELENGTH STEP SIZE FOR" in str(e):
                    logging.error(str(e))
                    logging.warn("Re-executing MOOG without the problematic line (it is recommended to manually remove that line to improve execution time)")
                    wave_A = float(str(e).split()[-3])
                    bad_idx = np.argmin(np.abs(linemasks['wave_A'][~bad] - wave_A))
                    bad[np.where(~bad)[0][bad_idx]] = True
                else:
                    raise
            else:
                success = True
        if np.all(bad):
            raise Exception("MOOG: No abundances could be calculated")
        if np.any(bad):
            nlines = len(linemasks)
            spec_abund = __reconstruct(nlines, spec_abund, bad)
            absolute_abund = __reconstruct(nlines, absolute_abund, bad)
            x_over_h = __reconstruct(nlines, x_over_h, bad)
            x_over_fe = __reconstruct(nlines, x_over_fe, bad)
        return spec_abund, absolute_abund, x_over_h, x_over_fe
    elif code == "width":
        #return __width_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, tmp_dir=tmp_dir)
        success = False
        bad = linemasks['wave_A'] < 0 # All to false
        # WIDTH can fail for some lines and it stops the execution, automatically detect that
        # filter the line and re-execute
        while not success and np.any(~bad):
            try:
                if ignore is not None:
                    filtered_ignore = ignore[~bad]
                else:
                    filtered_ignore = None
                spec_abund, absolute_abund, x_over_h, x_over_fe =  __width_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks[~bad], isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=filtered_ignore, verbose=verbose, tmp_dir=tmp_dir)
            except Exception, e:
                # WIDTH ERROR: 515.8551 -2.379  3.5  108707.400  2.0*********************01
                if "WIDTH ERROR:" in str(e):
                    logging.warn("Re-executing WIDTH without the problematic line (it is recommended to manually remove that line to improve execution time)")
                    wave_nm = float(str(e).split()[2])
                    bad_idx = np.argmin(np.abs(linemasks['wave_nm'][~bad] - wave_nm))
                    bad[np.where(~bad)[0][bad_idx]] = True
                else:
                    raise
            else:
                success = True
        if np.all(bad):
            raise Exception("WIDTH: No abundances could be calculated")
        if np.any(bad):
            nlines = len(linemasks)
            spec_abund = __reconstruct(nlines, spec_abund, bad)
            absolute_abund = __reconstruct(nlines, absolute_abund, bad)
            x_over_h = __reconstruct(nlines, x_over_h, bad)
            x_over_fe = __reconstruct(nlines, x_over_fe, bad)
        return spec_abund, absolute_abund, x_over_h, x_over_fe
    elif code == "spectrum":
        return __spectrum_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, gui_queue=gui_queue, timeout=timeout, tmp_dir=tmp_dir)

def __reconstruct(nlines, spec_abund, bad):
    spec_abund_reconstructed = np.zeros(nlines)
    spec_abund_reconstructed[~bad] = spec_abund
    spec_abund_reconstructed[bad] = np.nan
    return spec_abund_reconstructed


def __spectrum_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, gui_queue=None, timeout=1800, tmp_dir=None):
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

    linemasks_file = __spectrum_write_abundance_lines(linemasks)
    atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, tmp_dir=tmp_dir)
    abundances_file = write_solar_abundances(abundances, tmp_dir=tmp_dir)
    num_measures = len(linemasks)
    nlayers = len(atmosphere_layers)

    if ignore is None:
        ignore = np.ones(num_measures) # Compute fluxes for all the wavelengths
    # Filter out lines not supported by turbospectrum:
    lcode = linemasks['spectrum_support'] == "T"
    ignore[~lcode] = 0 # Ignore also lines not supported by spectrum

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__spectrum_determine_abundances_internal, args=(process_communication_queue, atmosphere_layers_file, linemasks_file, num_measures, ignore, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'verbose': verbose})
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

                # If the reference solar abundances given to SPECTRUM were
                # enhanced, we should correct it to put in the correct solar scale
                #abundances = enhance_solar_abundances(abundances, -1.*alpha)
                x_over_h = __correct_enhance_solar_abundances(linemasks, x_over_h, alpha)

                #import pudb
                #pudb.set_trace()
                #sun_log_Nx_over_Ntotal = abundances['Abund'][abundances['code'] == 26]]
                #x_absolute = free_abundances['Abund'][i] + 12. - sun_log_Nh_over_Ntotal # absolute, A(X)
                ##x_over_fe = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
                ##x_over_h = x_over_fe + self.MH()
                #x_over_h = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
                #x_over_fe = x_over_h - self.MH()
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

def __spectrum_determine_abundances_internal(process_communication_queue, atmosphere_model_file, linelist_file, num_measures, ignore, abundances_file, microturbulence_vel = 2.0, nlayers=56, verbose=0):
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
    if not is_spectrum_support_enabled():
        raise Exception("SPECTRUM support is not enabled")

    import synthesizer
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    abundances = synthesizer.abundances(atmosphere_model_file, linelist_file, num_measures, ignore, abundances_file, microturbulence_vel, nlayers, verbose, update_progress_func)
    process_communication_queue.put(abundances)



def __turbospectrum_read_abund_results(abundances_filename):
    """
    Read a turbospectrum resulting abundances file
    """
    # Fe  1      4802.880  3.64 -1.514   60.900   60.900 +-  1.000   0.60    8.023    8.054    8.085 strong? wings cut?
    f = "(?:[+-]?\d+\.\d+|NaN|[+-]?Inf)" # Float
    #
    atomic_line_pattern = re.compile("^\s*(\w+)\s+(\d+)\s+("+f+")\s+("+f+")\s+("+f+")\s+((?:"+f+"|\*+))\s+("+f+")\s+\+-\s+("+f+")\s+("+f+")\s+("+f+")\s+("+f+")\s+("+f+")(.*)\n$")
    nlines_read = 0
    data = []
    with open(abundances_filename, "ro") as f:
        for line_number, line in enumerate(f):
            if nlines_read >= 2:
                atomic_line = atomic_line_pattern.match(line)
                if atomic_line:
                    element, ion, wave_a, lower_state_eV, loggf, ew0, ew, ew_err, x_over_h, lower_abund, absolute_abund, upper_abund, comment = atomic_line.groups()
                    element = element + " " + ion
                    wave_a, lower_state_eV, loggf, ew, ew_err, \
                            x_over_h, lower_abund, absolute_abund, upper_abund = map(float, (wave_a, lower_state_eV, \
                            loggf, ew, ew_err, x_over_h, lower_abund, absolute_abund, upper_abund))
                    data.append((element, wave_a, wave_a/10., lower_state_eV, loggf, ew, ew_err, \
                            x_over_h, lower_abund, absolute_abund, upper_abund, comment))
                else:
                    raise Exception("Wrong format (line %i): %s" % (line_number, line))
            nlines_read += 1

    if len(data) == 0:
        abundances = data
    else:
        abundances = np.rec.fromarrays(zip(*data), dtype=[ \
                ('element', '|S4'), \
                ('wave_A', '<f8'), ('wave_nm', '<f8'), \
                ('lower_state_eV', float), \
                ('loggf', '<f8'), \
                ('ew', float), \
                ('ew_err', float), \
                ('x_over_h', float), \
                ('lower absolute abundance', float), \
                ('absolute abundance', float), \
                ('upper absolute abundance', float), \
                ('comment', '|S100') \
                ])
    return abundances


def __turbospectrum_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, tmp_dir=None):
    if not is_turbospectrum_support_enabled():
        raise Exception("Turbospectrum support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    turbospectrum_dir = ispec_dir + "/synthesizer/turbospectrum/"
    turbospectrum_data = turbospectrum_dir + "/DATA/"
    turbospectrum_eqwidt_lu = turbospectrum_dir + "bin/eqwidt_lu"

    radius = atmosphere_layers[0][-1]
    if radius > 2.0: # Compare to 2.0 instead of 1.0 to avoid floating point imprecisions
        spherical_model = True
    else:
        spherical_model = False

    atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, code="turbospectrum", atmosphere_filename=None, tmp_dir=tmp_dir)

    linemasks = linemasks.copy()
    # Filter out lines not supported by turbospectrum:
    lcode = linemasks['turbospectrum_support'] == "T"
    if ignore is not None:
        # Add specific lines to be ignored
        ignore = np.logical_or(np.logical_not(lcode), ignore == 0)
    else:
        ignore = np.logical_not(lcode)
    # Turbospectrum does not calculate abundances for lines with equivalent width equal to zero
    linemasks['ew'][ignore] = 0
    linemasks['ewr'][ignore] = 0

    wave_base = np.min(linemasks['wave_nm'])
    wave_top = np.max(linemasks['wave_nm'])
    wave_step = 0.001

    # Remove alpha enhancement to obtain the original solar scale
    original_abundances = enhance_solar_abundances(abundances, -1.*alpha)

    # Turbospectrum is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    atom_abundances = abundances[abundances['code'] <= 92]
    if len(atom_abundances) != 92:
        raise Exception("No abundances for all 92 elements!")
    efilter = np.logical_and(atom_abundances['code'] != 1, atom_abundances['code'] != 2)
    atom_abundances['Abund'][efilter] += MH

    opacities_filename = calculate_opacities(atmosphere_layers_file, atom_abundances, MH, microturbulence_vel, wave_base, wave_top, wave_step, verbose=verbose, opacities_filename=None)

    # For the determination of abundances, turbospectrum uses the atomic abundances
    # only to calculate [X/H], thus we have to provide the original non-modified
    # solar abundances
    atom_abundances = original_abundances[abundances['code'] <= 92]

    idx = []
    absolute_abund = []
    x_over_h = []
    for element in np.unique(linemasks['element']):
        efilter = linemasks['element'] == element
        idx = np.hstack((idx, np.where(efilter)[0])) # To recover later the order
        sublinemasks = linemasks[efilter]
        linelist_filename = write_atomic_linelist(sublinemasks, linelist_filename=None, code="turbospectrum", tmp_dir=tmp_dir)

        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False, dir=tmp_dir)
        out.close()
        abundances_filename = out.name

        # Temporary dir
        tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
        os.symlink(turbospectrum_data, tmp_execution_dir+"/DATA")
        previous_cwd = os.getcwd()
        os.chdir(tmp_execution_dir)

        command = turbospectrum_eqwidt_lu
        command_input = "'LAMBDA_MIN:'  '"+str(wave_base*10.)+"'\n"
        command_input += "'LAMBDA_MAX:'  '"+str(wave_top*10.)+"'\n"
        command_input += "'LAMBDA_STEP:' '"+str(wave_step*10.)+"'\n"
        command_input += "'INTENSITY/FLUX:' 'Flux'\n"
        command_input += "'COS(THETA)    :' '1.00'\n"
        command_input += "'ABFIND        :' '.true.'\n"
        command_input += "'MODELOPAC:' '"+opacities_filename+"'\n"
        command_input += "'RESULTFILE :' '"+abundances_filename+"'\n"
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
        if isotopes is not None:
            command_input += "'ISOTOPES : ' '"+str(len(isotopes))+"'\n"
            for isotope in isotopes:
                command_input += "%i.%03i  %.3f\n" % (isotope['atomic_code'], isotope['mass_number'], isotope['relative_abundance_in_the_solar_system'])
        command_input += "'NFILES   :' '1'\n"
        command_input += linelist_filename + "\n"
        if spherical_model:
            command_input += "'SPHERICAL:'  'T'\n"
        else:
            command_input += "'SPHERICAL:'  'F'\n"
        command_input += "  30\n"
        command_input += "  300.00\n"
        command_input += "  15\n"
        command_input += "  1.30\n"

        if verbose == 1:
            proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
        else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for the process to terminate
        out, err = proc.communicate(input=command_input)
        errcode = proc.returncode

        os.chdir(previous_cwd)

        resulting_abundances = __turbospectrum_read_abund_results(abundances_filename)

        # Turbospectrum does not return abundance for some elements sometimes, we should search for them
        # and add NaN for coherence
        absolute_abund_tmp = []
        x_over_h_tmp = []
        i = 0
        j = 0
        while i < len(sublinemasks) and j < len(resulting_abundances):
            equal_ew = np.abs(sublinemasks['ew'][i] - resulting_abundances['ew'][j]) < 0.1
            equal_wave_nm = np.abs(sublinemasks['wave_nm'][i] - resulting_abundances['wave_nm'][j]) < 0.002
            if equal_ew and equal_wave_nm:
                absolute_abund_tmp.append(resulting_abundances['absolute abundance'][j])
                x_over_h_tmp.append(resulting_abundances['x_over_h'][j])
                i += 1
                j += 1
            else:
                absolute_abund_tmp.append(np.nan)
                x_over_h_tmp.append(np.nan)
                i += 1
        while i < len(sublinemasks):
            absolute_abund_tmp.append(np.nan)
            x_over_h_tmp.append(np.nan)
            i += 1

        absolute_abund = np.hstack((absolute_abund, np.asarray(absolute_abund_tmp)))
        x_over_h = np.hstack((x_over_h, np.asarray(x_over_h_tmp)))

        os.remove(abundances_filename)
        os.remove(linelist_filename)
        shutil.rmtree(tmp_execution_dir)
    os.remove(atmosphere_layers_file)

    # Return abundances with the same order that the input linemasks
    sorted_idx_idx = np.argsort(idx)
    absolute_abund = absolute_abund[sorted_idx_idx]
    x_over_h = x_over_h[sorted_idx_idx]
    spec_abund = absolute_abund - 12.036 # For compatibility with SPECTRUM
    x_over_fe = x_over_h - MH

    return spec_abund, absolute_abund, x_over_h, x_over_fe


def __moog_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, tmp_dir=None):
    if not is_moog_support_enabled():
        raise Exception("MOOG support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    moog_dir = ispec_dir + "/synthesizer/moog/"
    moog_executable = moog_dir + "MOOGSILENT"

    tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
    #os.symlink(moog_dir, tmp_execution_dir+"/DATA")
    os.makedirs(tmp_execution_dir+"/DATA/")
    atmosphere_filename = tmp_execution_dir + "/model.in"
    linelist_file = tmp_execution_dir + "/lines.in"
    barklem_linelist_file = tmp_execution_dir + "/DATA/Barklem.dat"

    atmosphere_filename = write_atmosphere(atmosphere_layers, teff, logg, MH, code="moog", atmosphere_filename=atmosphere_filename, tmp_dir=tmp_dir)

    sorted_idx = np.argsort(linemasks, order='spectrum_moog_species')
    tmp_linemasks = linemasks.copy()

    # Filter out lines not supported by turbospectrum:
    lcode = linemasks['moog_support'] == "T"
    if ignore is not None:
        # Add specific lines to be ignored
        ignore = np.logical_or(np.logical_not(lcode), ignore == 0)
    else:
        ignore = np.logical_not(lcode)
    # Set ew to zero for lines to be ignored
    tmp_linemasks['ew'][ignore] = 0
    tmp_linemasks['ewr'][ignore] = 0

    sorted_tmp_linemasks = tmp_linemasks[sorted_idx]
    # MOOG does not support zero equivalent widths, so we will filter them:
    filtered = np.logical_or(sorted_tmp_linemasks['ew'] < 1e-9, np.isnan(sorted_tmp_linemasks['ew']))

    # Write (MOOG requires a separate file for damping coeff if we want to provide rad coeff and alpha from ABO theory)
    linelist_filename = write_atomic_linelist(sorted_tmp_linemasks[~filtered], linelist_filename=linelist_file, code="moog", tmp_dir=tmp_dir)
    linelist_filename = write_atomic_linelist(sorted_tmp_linemasks[~filtered], linelist_filename=barklem_linelist_file, code="moog_barklem", tmp_dir=tmp_dir)

    # Remove alpha enhancement to obtain the original solar scale
    original_abundances = enhance_solar_abundances(abundances, -1.*alpha)

    # MOOG is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    abundances = abundances.copy()
    efilter = np.logical_and(abundances['code'] != 1, abundances['code'] != 2)
    efilter = np.logical_and(efilter, abundances['code'] <= 92)
    abundances['Abund'][efilter] += MH

    # Append microturbulence, solar abundances and metallicity
    moog_atmosphere = open(atmosphere_filename, "a")
    atom_abundances = abundances[efilter] # Don't update hydrogen or helium abundances
    moog_atmosphere.write("  %.2f\n" % (microturbulence_vel))
    moog_atmosphere.write("NATOMS=   %i %.2f\n" % (len(atom_abundances), MH))
    for atom_abundance in atom_abundances:
        abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
        moog_atmosphere.write("%i  %.2f\n" % (atom_abundance['code'], abund))
    # Molecule list as used by Jorge Melendez (private communication)
    moog_atmosphere.write("NMOL      28\n")
    moog_atmosphere.write("  101.0   106.0   107.0   108.0   112.0  126.0\n")
    moog_atmosphere.write("  606.0   607.0   608.0\n")
    moog_atmosphere.write("  707.0   708.0\n")
    moog_atmosphere.write("  808.0   812.0   822.0   823.0   840.0\n")
    moog_atmosphere.write("  10108.0 10820.0 60808.0\n")
    moog_atmosphere.write("  6.1     7.1     8.1   12.1  20.1  22.1  23.1  26.1  40.1\n")
    moog_atmosphere.close()

    par_file = open(tmp_execution_dir + "/batch.par", "w")
    par_file.write("abfind\n")
    par_file.write("standard_out moog.std\n")
    par_file.write("summary_out  moog.sum\n")
    par_file.write("model_in     model.in\n")
    par_file.write("lines_in     lines.in\n")
    par_file.write("atmosphere   1\n") # controls the output of atmosphere quantities (print out the standard things about an atmsophere)
    par_file.write("molecules    1\n") # controls the molecular equilibrium calculations (1 = do molecular equilibrium but do not print results)
    par_file.write("units        0\n") # controls the units in which moog outputs the final spectrum (0 = angs, 1 = microns, 2 = 1/cm)
    par_file.write("lines        1\n") # controls the output of line data (print out standard information about the input line list)
    par_file.write("flux/int     0\n") # choses integrated flux or central intensity (0 = integrated flux calculations)
    par_file.write("damping      1\n") # Use Barklem.dat file with waals, alpha and rad damping coeff, if not found then do like damping = 0
    #par_file.write("damping      0\n") # Use single gamma damping coefficient (van der Waals) if provided, if zero then use Unsold equation
    par_file.write("freeform     0\n")  # Linelist format of 7 columns with numbers %10.3f and comment %10s
    par_file.write("plot         0\n")
    par_file.close()

    previous_cwd = os.getcwd()
    os.chdir(tmp_execution_dir)

    command = moog_executable
    command_input = ""

    # MOOG clears the screen, so it is better to always don't do verbose this part
    # also the information printed is not specially useful
    #if verbose == 1:
        #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
    #else:
    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait for the process to terminate
    out, err = proc.communicate(input=command_input)
    errcode = proc.returncode

    for line in out.split("\n"):
        if "I QUIT" in line:
            # MOOG ERROR: CANNOT DECIDE ON LINE WAVELENGTH STEP SIZE FOR   5158.62   I QUIT!
            logging.error("MOOG ERROR: %s" % (line))
            os.chdir(previous_cwd)
            shutil.rmtree(tmp_execution_dir)
            raise Exception("MOOG ERROR: %s" % (line))

    absolute_abund = []
    x_over_h = []
    results = open(tmp_execution_dir + "/moog.sum", "r")
    results_lines = results.readlines()
    # Abundance Results for Species Fe I         (input abundance =   7.450)
    #   4893.813   26.10000   2.828  -4.267    23.62    -5.316     7.650    0.018
    f = "(?:[+-]?\d+\.\d+|NaN|[+-]?Inf)" # Float
    atomic_line_pattern = re.compile("^\s+("+f+")\s+(.*)\n$")
    line_number = 0
    for i, line in enumerate(results_lines):
        atomic_line = atomic_line_pattern.match(line)
        if atomic_line:
            values = map(float, line.split())
            # Make sure we recover the abundance for the good line by checking the wavelength
            wave = values[0] # A
            while filtered[line_number] or np.abs(sorted_tmp_linemasks['wave_A'][line_number] - wave) > 0.001:
                absolute_abund.append(np.nan)
                x_over_h.append(np.nan)
                if filtered[line_number]:
                    #logging.warn("Missed line %.3f because of NaN equivalent width, not supported by MOOG or requested to be ignored" % (sorted_tmp_linemasks['wave_A'][line_number]))
                    pass
                else:
                    logging.warn("Missed line %.3f" % (sorted_tmp_linemasks['wave_A'][line_number], wave))
                line_number += 1
            absolute_abund.append(values[-2])
            #x_over_h.append(values[-1]) # This MOOG value is not consistent with absolute_abund - solar abundance
            if reference_abund is None:
                species = int(float(values[1])) # Convert from '26.0' or '26.1' to 26
                reference_abund = original_abundances['Abund'][original_abundances['code'] == species][0] + 12.036
            x_over_h.append(values[-2] - reference_abund)
            line_number += 1
        else:
            # Detect new species results
            if "Abundance Results for Species" in line:
                reference_abund = None
    for last_missing_lines in sorted_tmp_linemasks['wave_A'][line_number:]:
        absolute_abund.append(np.nan)
        x_over_h.append(np.nan)

    absolute_abund = np.asarray(absolute_abund)
    x_over_h = np.asarray(x_over_h)

    os.chdir(previous_cwd)
    shutil.rmtree(tmp_execution_dir)

    # Return abundances with the same order that the input linemasks
    sorted_idx_idx = np.argsort(sorted_idx)
    absolute_abund = absolute_abund[sorted_idx_idx]
    x_over_h = x_over_h[sorted_idx_idx]
    spec_abund = absolute_abund - 12.036 # For compatibility with SPECTRUM
    x_over_fe = x_over_h - MH

    return spec_abund, absolute_abund, x_over_h, x_over_fe


def __width_determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, isotopes, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, tmp_dir=None):
    if not is_width_support_enabled():
        raise Exception("WIDTH support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    atmos_dir = ispec_dir + "/synthesizer/atmos/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        width_executable = atmos_dir + "bin.amd64/width9.exe"
    else:
        width_executable = atmos_dir + "bin.ia32/width9.exe"
    width_molecules = atmos_dir + "lines/molecules.dat"

    tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
    os.symlink(width_molecules, tmp_execution_dir+"/fort.2")

    # Remove alpha enhancement to obtain the original solar scale
    original_abundances = enhance_solar_abundances(abundances, -1.*alpha)

    previous_cwd = os.getcwd()
    os.chdir(tmp_execution_dir)

    linemasks = linemasks.copy()
    # Filter out lines not supported by width:
    lcode = linemasks['width_support'] == "T"
    if ignore is not None:
        # Add specific lines to be ignored
        ignore = np.logical_or(np.logical_not(lcode), ignore == 0)
    else:
        ignore = np.logical_not(lcode)
    # Zero EW will be filtered out later
    linemasks['ew'][ignore] = 0
    linemasks['ewr'][ignore] = 0

    command = width_executable
    command_input = ""
    command_input += "VTUR\n"
    command_input += "    1 %.2f\n" % (microturbulence_vel)
    lines_input = ""
    filtered = []
    for line in linemasks:
        ew_picometer = line['ew'] / 10.
        if line['width_support'] == "False" or np.isnan(line['ew']) or line['ew'] < 1e-9:
            filtered.append(True)
            continue
        else:
            lower_state = line['lower_state_cm1']
            upper_state = line['upper_state_cm1']
            lower_level = line['lower_j']
            upper_level = line['upper_j']
            width_species = line['width_species']
            nelion = 0 # ?
            waals = line['waals_single_gamma_format']
            alpha = 0.00 # Omara theory (ABO) sigma.alpha (sigma was transformed to waals_single_gamma_format) and alpha is provided separately
            if line['spectrum_transition_type'] == "AO":
                alpha = line['waals'] % 1 # Decimal part
            if line['rad'] == 0 and line['stark'] == 0  and waals == 0:
                # WIDTH does not all damping parameters to zero
                filtered.append(True)
            else:
                filtered.append(False)
                command_input += "LINE     %.2f  %.4f    STARNAME\n" % (ew_picometer, line['mu'])
                command_input += "  %.4f %.3f  %.1f   %.3f  %.1f  %.3f     %s\n" % (line['wave_nm'], line['loggf'], lower_level, lower_state, upper_level, upper_state, width_species)
                command_input += "  %.4f  %i  %.2f %.2f  %.2f     0  0  0.000  0  0.000    0    0   0  0 %.3f\n" % (line['wave_nm'], nelion, line['rad'], line['stark'], waals, alpha)

    command_input += "END\n"

    command_input += "TEFF   %.0f  GRAVITY %.5f LTE\n" % (teff, logg)
    # command_input += "TITLE  [0.0] VTURB=1.0 KM/SEC  L/H=1.25 NOVER NEW ODF ASPLUND ABUNDANCES\n"
    # command_input += " OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0\n"
    # mixing_length_param = 1.25
    # command_input += " CONVECTION ON   %.2f TURBULENCE OFF  0.00  0.00  0.00  0.00\n" % (mixing_length_param)
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
    command_input += "END\n"
    command_input += "STOP\n"

    # Never verbose because WIDTH's output is printed on stdout (not saved on a file)
    # and it is needed by iSpec
    #if verbose == 1:
        #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
    #else:
    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    # wait for the process to terminate
    out, err = proc.communicate(input=command_input)
    errcode = proc.returncode
    #c = open("conf.txt", "w")
    #c.writelines(command_input)
    #c.close()
    #import pudb
    #pudb.set_trace()

    line_number = 0
    not_converged = False
    spec_abund = np.zeros(len(linemasks))
    absolute_abund = np.zeros(len(linemasks))
    x_over_h = np.zeros(len(linemasks))
    x_over_fe = np.zeros(len(linemasks))
    # Default to nan (needed for cases where the last lines are not computed/reported by WIDTH):
    spec_abund[:] = np.nan
    absolute_abund[:] = np.nan
    x_over_h[:] = np.nan
    x_over_fe[:] = np.nan
    #  666.7710 -2.112  1.0   36972.005  0.0   51966.006    26.00
    #  666.7710   0  8.09 -5.63 -7.67     0 0  0  0.000  0  0.000  0    0
    #     STARNAME                                                                  666.7724  8.09 -5.63 -7.67  0.60      1.21    -4.391
    #         VTURB  ABUND    -4.59    -4.38    -4.09    -3.59    -4.39
    #          1.00     EW   -0.090    0.093    0.320    0.624    0.083
    #                          0.81     1.24     2.09     4.21     1.21
    #                DEPTH     0.61     0.59     0.55     0.40     0.60
    #
    #  667.7985 -1.418  1.0   21712.003  0.0   36690.005    26.00
    #  667.7985   0  8.07 -6.09 -7.64     0 0  0  0.000  0  0.000  0    0
    #0    STARNAME                                                               **********  8.07 -6.09 -7.64 -0.18     12.50    -4.729
    #         VTURB  ABUND    -5.59    -5.09    -4.74    -4.59    -4.73
    #          1.00     EW    0.848    0.993    1.095    1.140    1.097
    #                          7.04     9.84    12.43    13.80    12.50
    #                DEPTH     0.10    -0.09    -0.18    -0.20    -0.18
    out_lines = out.split("\n")
    for i, line in enumerate(out_lines):
        if "BETTER LUCK NEXT TIME" in line:
            logging.error("WIDTH" + out_lines[i-2])
            logging.error("WIDTH" + out_lines[i-1])
            logging.error("WIDTH" + out_lines[i])
            os.chdir(previous_cwd)
            shutil.rmtree(tmp_execution_dir)
            raise Exception("WIDTH ERROR: %s" % (out_lines[i-2]))
        if "NOT CONVERGED" in line:
            #  679.3258 -2.326  4.0   32875.157  4.0       5.900    26.00
            #  679.3258   0  7.56 -5.99 -7.82     0 0  0  0.000  0  0.000  0    0
            #
            #
            # NOT CONVERGED
            #     STARNAME                                                                  679.3262  7.56 -5.99 -7.82  0.60      1.33    -8.679
            #         VTURB  ABUND    -9.90    -8.68    -8.67    -8.58    -8.07    -6.59    -5.59    -4.59    -8.68
            #          1.00     EW   -1.029    0.125    0.133    0.208    0.569    1.055    1.337    1.761    0.124
            #                          0.09     1.33     1.36     1.61     3.71    11.35    21.71    57.65     1.33
            #                DEPTH     0.64     0.60     0.59     0.58     0.38    -0.49    -0.60    -0.20     0.60
            not_converged = True
        if "STARNAME" not in line:
            continue

        if not_converged:
            spec_abund[line_number] = np.nan
            absolute_abund[line_number] = np.nan
            x_over_h[line_number] = np.nan
            x_over_fe[line_number] = np.nan
            line_number += 1
            not_converged = False
            continue

        # Make sure we recover the abundance for the good line by checking the wavelength
        wave = float(out_lines[i-1].split()[0]) # nm
        while filtered[line_number] or np.abs(linemasks['wave_nm'][line_number] - wave) > 0.0001:
            spec_abund[line_number] = np.nan
            absolute_abund[line_number] = np.nan
            x_over_h[line_number] = np.nan
            x_over_fe[line_number] = np.nan
            if filtered[line_number]:
                #logging.warn("Missed line %.3f because of NaN equivalent width, bad atomic information or not supported by WIDTH" % (linemasks['wave_A'][line_number]))
                pass
            else:
                logging.warn("Missed line %.3f [%.3f]" % (linemasks['wave_A'][line_number], wave))
            line_number += 1

        values = line.split()
        spec_abund[line_number] = float(values[-1])
        absolute_abund[line_number] = spec_abund[line_number] + 12.036
        species = int(float(linemasks['width_species'][line_number])) # Convert from '26.0' or '26.1' to 26
        solar_abund = abundances['Abund'][abundances['code'] == species]
        x_over_h[line_number] = spec_abund[line_number] - solar_abund
        x_over_fe[line_number] = x_over_h[line_number] - MH
        line_number += 1

    os.chdir(previous_cwd)
    shutil.rmtree(tmp_execution_dir)

    return spec_abund, absolute_abund, x_over_h, x_over_fe


def __correct_enhance_solar_abundances(linemasks, abundances, alpha_enhancement):
    """
    Descales alpha elements and CNO abundances.
    """
    abundances = abundances.copy()
    code = map(int, map(float, linemasks['spectrum_moog_species'])) # Convert from '26.0' or '26.1' to 26

    #  6|C|Carbon|14|2|6|4|4
    c = code == 6
    #  7|N|Nitrogen|15|2|7|6|8
    n = code == 7
    #  8|O|Oxygen|16|2|8|3|3
    o = code == 8

    # 10|Ne|Neon|18|2|10|5|5
    # 12|Mg|Magnesium|2|3|12|7|6
    alpha = np.logical_or(code == 10, code == 12)
    # 14|Si|Silicon|14|3|14|8|7
    alpha = np.logical_or(alpha, code == 14)
    # 16|S|Sulfur|16|3|16|10|9
    alpha = np.logical_or(alpha, code == 16)
    # 18|Ar|Argon|18|3|18|14|11
    alpha = np.logical_or(alpha, code == 18)
    # 20|Ca|Calcium|2|4|20|12|12
    alpha = np.logical_or(alpha, code == 20)
    # 22|Ti|Titanium|4|4|22|22|19
    alpha = np.logical_or(alpha, code == 22)

    #abundances[c] -= c_enhancement
    #abundances[n] -= n_enhancement
    #abundances[o] -= o_enhancement
    abundances[o] -= alpha_enhancement # MARCS standard enhances Oxygen as well as alpha elements at the same rate
    abundances[alpha] -= alpha_enhancement

    return abundances
