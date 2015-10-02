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
from lines import write_atomic_linelist
from common import is_turbospectrum_support_enabled, is_spectrum_support_enabled
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
from Queue import Empty
import subprocess
import shutil

import log
import logging
import re



def write_abundance_lines(linemasks, filename=None, tmp_dir=None):
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
    Calculates alpha_enhancement, c_enhancement, n_enhancement, o_enhancement
    from the metallicity (MH). By default the standard composition of MARCS
    model atmospheres is used, which is the most logical for most analysis.
    """
    if scale is None:
        # MARCS standard scale
        # - http://marcs.astro.uu.se/
        scale = [(+1.00,   0.00,    0.00,    0.00,    0.00,  ),
                 (+0.75,   0.00,    0.00,    0.00,    0.00,  ),
                 (+0.50,   0.00,    0.00,    0.00,    0.00,  ),
                 (+0.25,   0.00,    0.00,    0.00,    0.00,  ),
                 (+0.00,   0.00,    0.00,    0.00,    0.00,  ),
                 (-0.25,   +0.10,   0.00,    0.00,    +0.10, ),
                 (-0.50,   +0.20,   0.00,    0.00,    +0.20, ),
                 (-0.75,   +0.30,   0.00,    0.00,    +0.30, ),
                 (-1.00,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-1.50,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-2.00,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-2.50,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-3.00,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-4.00,   +0.40,   0.00,    0.00,    +0.40, ),
                 (-5.00,   +0.40,   0.00,    0.00,    +0.40, )]
        scale = np.array(scale, dtype=[('[Fe/H]', float), ('[a/Fe]', float), ('[C/Fe]', float), ('[N/Fe]', float), ('[O/Fe]', float)])
    scale.sort(order='[Fe/H]') # It should be sorted or interpolation does not work properly
    alpha_enhancement = np.interp(MH, scale['[Fe/H]'], scale['[a/Fe]'])
    c_enhancement = np.interp(MH, scale['[Fe/H]'], scale['[C/Fe]'])
    n_enhancement = np.interp(MH, scale['[Fe/H]'], scale['[N/Fe]'])
    o_enhancement = np.interp(MH, scale['[Fe/H]'], scale['[O/Fe]'])
    return alpha_enhancement, c_enhancement, n_enhancement, o_enhancement




def enhance_solar_abundances(abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement):
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

    abundances['Abund'][c] += c_enhancement
    abundances['Abund'][n] += n_enhancement
    abundances['Abund'][o] += o_enhancement
    abundances['Abund'][alpha] += alpha_enhancement

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


def determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, gui_queue=None, timeout=1800, isotopes=None, code="spectrum", tmp_dir=None, enhance_abundances=True, scale=None):
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    if code == "turbospectrum":
        return __turbospectrum_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, tmp_dir=tmp_dir, enhance_abundances=enhance_abundances, scale=scale)
    elif code == "moog":
        return __moog_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, isotopes, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, tmp_dir=tmp_dir, enhance_abundances=enhance_abundances, scale=scale)
    elif code == "spectrum":
        return __spectrum_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = microturbulence_vel, ignore=ignore, verbose=verbose, gui_queue=gui_queue, timeout=timeout, tmp_dir=tmp_dir, enhance_abundances=enhance_abundances, scale=scale)

def __spectrum_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, gui_queue=None, timeout=1800, tmp_dir=None, enhance_abundances=True, scale=None):
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

    # Enhance alpha elements + CNO abundances following MARCS standard composition
    if enhance_abundances:
        alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(MH, scale=scale)
        abundances = enhance_solar_abundances(abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

    linemasks_file = write_abundance_lines(linemasks)
    atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH, tmp_dir=tmp_dir)
    abundances_file = write_solar_abundances(abundances, tmp_dir=tmp_dir)
    num_measures = len(linemasks)
    nlayers = len(atmosphere_layers)
    if ignore == None:
        ignore = np.ones(num_measures) # Compute fluxes for all the wavelengths

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__determine_abundances, args=(process_communication_queue, atmosphere_layers_file, linemasks_file, num_measures, ignore, abundances_file,), kwargs={'microturbulence_vel': microturbulence_vel, 'nlayers': nlayers, 'verbose': verbose})
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
                if enhance_abundances:
                    x_over_h = __correct_enhance_solar_abundances(linemasks, x_over_h, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)
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

def __determine_abundances(process_communication_queue, atmosphere_model_file, linelist_file, num_measures, ignore, abundances_file, microturbulence_vel = 2.0, nlayers=56, verbose=0):
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
    abundances = np.rec.fromarrays(zip(*data), dtype=[ \
            ('element', '|S4'), \
            ('wave (A)', '<f8'), ('wave (nm)', '<f8'), \
            ('lower state (eV)', float), \
            ('log(gf)', '<f8'), \
            ('ew', float), \
            ('ew_err', float), \
            ('x_over_h', float), \
            ('lower absolute abundance', float), \
            ('absolute abundance', float), \
            ('upper absolute abundance', float), \
            ('comment', '|S100') \
            ])
    return abundances


def __turbospectrum_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, isotopes, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, tmp_dir=None, enhance_abundances=True, scale=None):
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
    if ignore is not None:
        # Turbospectrum does not calculate abundances for lines with equivalent width equal to zero
        linemasks = linemasks.copy()
        linemasks['ew'][ignore == 0] = 0
        linemasks['ewr'][ignore == 0] = 0

    wave_base = np.min(linemasks['wave (nm)'])
    wave_top = np.max(linemasks['wave (nm)'])
    wave_step = 0.001

    # Enhance alpha elements + CNO abundances following MARCS standard composition
    original_abundances = abundances.copy()
    if enhance_abundances:
        alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(MH, scale=scale)
        abundances = enhance_solar_abundances(abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

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
            equal_wave_nm = np.abs(sublinemasks['wave (nm)'][i] - resulting_abundances['wave (nm)'][j]) < 0.002
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


def __moog_determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, isotopes, abundances, microturbulence_vel = 2.0, ignore=None, verbose=0, tmp_dir=None, enhance_abundances=True, scale=None):
    if not is_turbospectrum_support_enabled():
        raise Exception("Turbospectrum support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../"
    moog_dir = ispec_dir + "/synthesizer/moog/"
    moog_executable = moog_dir + "MOOGSILENT"

    tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
    os.symlink(moog_dir, tmp_execution_dir+"/DATA")
    atmosphere_filename = tmp_execution_dir + "/model.in"
    linelist_file = tmp_execution_dir + "/lines.in"

    sorted_idx = np.argsort(linemasks, order='species')
    atmosphere_filename = write_atmosphere(atmosphere_layers, teff, logg, MH, code="moog", atmosphere_filename=atmosphere_filename, tmp_dir=tmp_dir)
    linelist_filename = write_atomic_linelist(linemasks[sorted_idx], linelist_filename=linelist_file, code="moog", tmp_dir=tmp_dir)

    # Enhance alpha elements + CNO abundances following MARCS standard composition
    original_abundances = abundances.copy()
    if enhance_abundances:
        alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = determine_abundance_enchancements(MH, scale=scale)
        abundances = enhance_solar_abundances(abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

    # Append microturbulence, solar abundances and metallicity
    moog_atmosphere = open(atmosphere_filename, "a")
    atom_abundances = abundances[abundances['code'] <= 92]
    moog_atmosphere.write("  %.2f\n" % (microturbulence_vel))
    moog_atmosphere.write("NATOMS=   %i %.2f\n" % (len(atom_abundances), MH))
    for atom_abundance in atom_abundances:
        abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
        moog_atmosphere.write("%i  %.2f\n" % (atom_abundance['code'], abund))
    moog_atmosphere.write("NMOL      22\n")
    moog_atmosphere.write("  101.0   106.0   107.0   108.0   112.0  126.0\n")
    moog_atmosphere.write("  606.0   607.0   608.0\n")
    moog_atmosphere.write("  707.0   708.0\n")
    moog_atmosphere.write("  808.0   812.0   822.0\n")
    moog_atmosphere.write("  10108.0 60808.0\n")
    moog_atmosphere.write("  6.1     7.1     8.1   12.1  22.1  26.1")
    moog_atmosphere.close()

    par_file = open(tmp_execution_dir + "/batch.par", "w")
    par_file.write("abfind\n")
    par_file.write("standard_out moog.std\n")
    par_file.write("summary_out  moog.sum\n")
    par_file.write("model_in     model.in\n")
    par_file.write("lines_in     lines.in\n")
    par_file.write("atmosphere   1\n")
    par_file.write("molecules    1\n")
    par_file.write("lines        1\n")
    par_file.write("flux/int     0\n")
    par_file.write("damping      1\n")
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


    absolute_abund = []
    x_over_h = []
    results = open(tmp_execution_dir + "/moog.sum", "r")
    results_lines = results.readlines()
    # Abundance Results for Species Fe I         (input abundance =   7.450)
    #   4893.813   26.10000   2.828  -4.267    23.62    -5.316     7.650    0.018
    f = "(?:[+-]?\d+\.\d+|NaN|[+-]?Inf)" # Float
    atomic_line_pattern = re.compile("^\s+("+f+")\s+(.*)\n$")
    for line in results_lines:
        atomic_line = atomic_line_pattern.match(line)
        if atomic_line:
            values = map(float, line.split())
            absolute_abund.append(values[-2])
            #x_over_h.append(values[-1]) # This MOOG value is not consistent with absolute_abund - solar abundance
            if reference_abund is None:
                species = int(float(values[1])) # Convert from '26.0' or '26.1' to 26
                reference_abund = original_abundances['Abund'][original_abundances['code'] == species][0] + 12.036
            x_over_h.append(values[-2] - reference_abund)
        else:
            # Detect new species results
            if "Abundance Results for Species" in line:
                reference_abund = None

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


def __correct_enhance_solar_abundances(linemasks, abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement):
    """
    Descales alpha elements and CNO abundances.
    """
    abundances = abundances.copy()
    code = map(int, map(float, linemasks['species'])) # Convert from '26.0' or '26.1' to 26

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

    abundances[c] -= c_enhancement
    abundances[n] -= n_enhancement
    abundances[o] -= o_enhancement
    abundances[alpha] -= alpha_enhancement

    return abundances
