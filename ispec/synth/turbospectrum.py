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
import sys
import numpy as np
import subprocess
import shutil
import re
import glob
import tempfile
import logging

from ispec.atmospheres import write_atmosphere, calculate_opacities
from ispec.lines import write_atomic_linelist
from ispec.common import which, is_turbospectrum_support_enabled
from ispec.spectrum import create_spectrum_structure, resample_spectrum
from .effects import _filter_linelist, apply_post_fundamental_effects


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, use_molecules=False, nlte_departure_coefficients=None, tmp_dir=None, timeout=1800):
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, use_molecules=use_molecules, nlte_departure_coefficients=nlte_departure_coefficients, tmp_dir=tmp_dir, timeout=timeout)

def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, use_molecules=False, nlte_departure_coefficients=None, tmp_dir=None, timeout=1800):
    if not is_turbospectrum_support_enabled():
        raise Exception("Turbospectrum support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
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
    linelist = _filter_linelist(linelist, regions)


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

    is_marcs_model = len(atmosphere_layers[0]) == 11
    opacities_file = calculate_opacities(atmosphere_layers_file, atom_abundances, MH, microturbulence_vel, global_wave_base-10, global_wave_top+10, wave_step, verbose=verbose, opacities_filename=None, tmp_dir=tmp_dir, is_marcs_model=is_marcs_model)

    if linelist_file is None:
        remove_tmp_linelist_file = True
        linelist_filename = write_atomic_linelist(linelist, linelist_filename=linelist_file, code="turbospectrum", tmp_dir=tmp_dir)
    else:
        linelist_filename = linelist_file

    #--------------------------------------------------------------------------------
    # NLTE
    #--------------------------------------------------------------------------------
    tmp_departures_dir = None
    absolute_atom_model_dirname = None
    nlte_info_file_content = []
    nlte_available = []
    nlte_ignored = []
    if nlte_departure_coefficients is not None:
        for element_name, (interpolated_departures, interpolated_tau, ndep, nk, atomic_number, absolute_abundance, absolute_atom_model_filename, ) in nlte_departure_coefficients.items():
            if atomic_number > 1: # hydrogen (atomic number == 1) is usually encoded into radiative transfer codes and not necessarily present in the atomic linelist
                # Check if there is any atomic line related to this element that has NLTE data
                element_lines_filter = np.floor(np.asarray(linelist['turbospectrum_species'], dtype=float)) == atomic_number
                element_linelist = linelist[element_lines_filter]
                linelist_has_nlte_data_for_element = 'T' in element_linelist['nlte'] and (np.any(element_linelist[element_linelist['nlte'] == 'T']['nlte_label_low'] != 'none') or np.any(element_linelist[element_linelist['nlte'] == 'T']['nlte_label_up'] != 'none'))
                if not linelist_has_nlte_data_for_element:
                    # Skip interpolattion of departure coefficient for elements that do not have any line with NLTE data in the linelist
                    nlte_ignored.append(element_name)
                    continue
            # interpolate
            if tmp_departures_dir is None:
                tmp_departures_dir = tempfile.mkdtemp(dir=tmp_dir)
            if absolute_atom_model_dirname is None:
                absolute_atom_model_dirname = os.path.dirname(absolute_atom_model_filename)
            else:
                assert absolute_atom_model_dirname == os.path.dirname(absolute_atom_model_filename)
            atom_model_filename = os.path.basename(absolute_atom_model_filename)
            out = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=tmp_departures_dir, encoding='utf-8')
            out.close()
            departure_coefficient_filename = out.name

            # Overwrite tau (depth) to avoid turbospectrum's error "tau scales differ in model atmos and in departure coefficient file"
            # - tau is interpolated at the atmosphere level and at the departure coefficient level
            #   both from the same MARCS model atmospheres, however, due to differences in how the
            #   interpolation is done, the interpolated tau can be slighly different (difference
            #   bigger than the hardcoded turbospectrum threshold)
            interpolated_tau = atmosphere_layers[:, 7]

            # Write departure coefficient file
            _turbospectrum_write_departure_coefficient_file(
                departure_coefficient_filename,
                ndep=ndep,
                nlevel=nk,
                tau_log10=interpolated_tau,
                depart_coeffs=interpolated_departures,
                atmos_str='interpolated',
                abundance=absolute_abundance,
                binary=False
            )
            nlte_available.append(element_name)
            nlte_info_file_content.append((atomic_number, element_name, atom_model_filename, os.path.basename(departure_coefficient_filename),))


    nlte_info_filename = None
    if len(nlte_info_file_content) > 0:
        out = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=tmp_departures_dir, encoding='utf-8')
        out.close()
        nlte_info_filename = out.name
        with open(nlte_info_filename, "w") as f:
            f.write("# This file controls which species are treated in LTE/NLTE\n")
            f.write("# It also gives the path to the model atom and the departure files\n")
            f.write("# First created 2021-02-22\n")
            f.write("# if a species is absent it is assumed to be LTE\n")
            f.write("#\n")
            f.write("# each line contains :\n")
            f.write("# atomic number / name / (n)lte / model atom / departure file / binary or ascii departure file\n")
            f.write("#\n")
            f.write("# path for model atom files     ! don't change this line !\n")
            f.write(f"{absolute_atom_model_dirname}/\n")
            f.write("#\n")
            f.write("# path for departure files      ! don't change this line !\n")
            f.write(f"{tmp_departures_dir}/\n")
            f.write("#\n")
            f.write("# atomic (N)LTE setup\n")
            for atomic_number, element_name, atom_model_filename, departure_coefficient_filename in nlte_info_file_content:
                f.write(f"{atomic_number}\t'{element_name}'\t'nlte'\t'{atom_model_filename}'\t'{departure_coefficient_filename}'\t'ascii'\n")
    #--------------------------------------------------------------------------------

    synth_fluxes = []
    synth_waveobs = []
    non_positive_result = False
    for segment in segments:
        wave_base = segment['wave_base']
        wave_top = segment['wave_top']
        # Temporary file
        out = tempfile.NamedTemporaryFile(mode="wt", delete=False, dir=tmp_dir, encoding='utf-8')
        out.close()
        synth_spectrum_filename = out.name

        # Temporary dir
        tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
        os.symlink(turbospectrum_data, tmp_execution_dir+"/DATA")
        os.symlink(molecules_dir, tmp_execution_dir+"/molecules")
        previous_cwd = os.getcwd()
        os.chdir(tmp_execution_dir)

        command = turbospectrum_bsyn_lu
        command_input = ""
        if nlte_info_filename is not None:
            command_input += "'NLTE :'          '.true.'\n"
            command_input += f"'NLTEINFOFILE:'  '{nlte_info_filename}'\n"
        else:
            command_input += "'NLTE :'          '.false.'\n"
        command_input += "'LAMBDA_MIN:'  '"+str(wave_base*10.)+"'\n"
        command_input += "'LAMBDA_MAX:'  '"+str(wave_top*10.)+"'\n"
        command_input += "'LAMBDA_STEP:' '"+str(wave_step*10.)+"'\n"
        if len(regions) > 1:
            segments_filename = "segments.txt"
            command_input += f"'SEGMENTSFILE:' '{segments_filename}'\n"
            with open(segments_filename, "w") as f:
                f.write("; synthesis segments\n")
                for region in regions:
                    f.write(f"{region['wave_base']*10.:13.2f}{region['wave_top']*10.:13.2f}\n")
        command_input += "'INTENSITY/FLUX:' 'Flux'\n"
        command_input += "'COS(THETA)    :' '1.00'\n"
        command_input += "'ABFIND        :' '.false.'\n"
        if is_marcs_model:
            command_input += "'MARCS-FILE    :' '.true.'\n"
        else:
            command_input += "'MARCS-FILE    :' '.false.'\n"
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
                name, file_wave_base, file_wave_top = re.match(r"(.*)_(\d+)-(\d+)\.bsyn", os.path.basename(filename)).groups()
                file_wave_base = float(file_wave_base)
                file_wave_top = float(file_wave_top)
                if (file_wave_base >= wave_base and file_wave_top <= wave_top) or \
                        (wave_base >= file_wave_base and wave_base <= file_wave_top ) or \
                        (wave_top >= file_wave_base and wave_top <= file_wave_top ):
                    molecules += filename + "\n"
                    num_molecules_files += 1
            command_input += "'NFILES   :' '%i'\n" % (2 + num_molecules_files)
            if num_molecules_files > 0:
                command_input += molecules
            else:
                print(f"No molecules found between '{wave_base}' and '{wave_top}'. Are there turbospectrum molecule files ('input/linelists/turbospectrum/molecules/*.bsyn')?")
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
        out, err = proc.communicate(input=command_input.encode('utf-8'))
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
            print(out)
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
    if tmp_departures_dir is not None:
        shutil.rmtree(tmp_departures_dir)

    synth_spectrum = create_spectrum_structure(synth_waveobs, synth_fluxes)
    synth_spectrum.sort(order=['waveobs'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(synth_spectrum['waveobs'], waveobs):
        synth_spectrum = resample_spectrum(synth_spectrum, waveobs, method="linear", zero_edges=True)

    segments = None
    vrad = (0,)
    synth_spectrum['flux'] = apply_post_fundamental_effects(synth_spectrum['waveobs'], synth_spectrum['flux'], segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)


    if verbose == 1:
        if nlte_info_filename is not None:
            if len(nlte_ignored) > 0:
                print(f"NLTE for elements: {nlte_available} | Not in region: {nlte_ignored}")
            else:
                print(f"NLTE for elements: {nlte_available}")
        else:
            print("Only LTE")

    return synth_spectrum['flux']

#--------------------------------------------------------------------------------
def _write_fortran_record(f, data):
    """
    Writes a Fortran-style unformatted record to a file.
    A record is data surrounded by 4-byte integers specifying the data length.
    """
    # '<i' means little-endian 4-byte integer
    len_bytes = struct.pack('<i', len(data))
    f.write(len_bytes)
    f.write(data)
    f.write(len_bytes)

def _turbospectrum_write_departure_coefficient_file(filename, ndep, nlevel, tau_log10, depart_coeffs, atmos_str='default_atmos', abundance=7.50, binary=False):
    """
    Creates a departure coefficient file in the format expected by the Fortran code.

    Args:
        filename (str): The name of the file to create.
        ndep (int): The number of depth points.
        nlevel (int): The number of energy levels.
        tau_log10 (np.ndarray): 1D array of log10(tau) values, shape (ndep,).
        depart_coeffs (np.ndarray): 2D array of departure coefficients.
                                   MUST have shape (nlevel, ndep).
        atmos_str (str): A descriptive string for the atmosphere model.
        abundance (float): The NLTE abundance to write to the file.
        binary (bool): If True, creates a binary file. If False, creates ASCII.
    """
    if depart_coeffs.shape != (nlevel, ndep):
        raise ValueError(f"Shape of depart_coeffs is {depart_coeffs.shape}, "
                         f"but expected ({nlevel}, {ndep})")

    if binary:
        # Create a little-endian, unformatted binary file
        #print(f"Creating BINARY departure file: {filename}")
        with open(filename, 'wb') as f:
            # Record 1: header_dep1 (500 bytes)
            header1 = atmos_str.ljust(500)
            _write_fortran_record(f, header1.encode('latin-1'))

            # Record 2: abundance_nlte (4-byte float)
            _write_fortran_record(f, struct.pack('<f', abundance))

            # Record 3: header_dep2 (1000 bytes)
            header2 = "Binary departure file created with Python".ljust(1000)
            _write_fortran_record(f, header2.encode('latin-1'))

            # Record 4: ndepth_read (4-byte integer)
            _write_fortran_record(f, struct.pack('<i', ndep))

            # Record 5: modnlevel_read (4-byte integer)
            _write_fortran_record(f, struct.pack('<i', nlevel))

            # Record 6: taumod array (linear tau)
            tau_linear = 10.0**tau_log10
            # '<{ndep}f' creates a format string like '<56f' for 56 floats
            tau_data = struct.pack(f'<{ndep}f', *tau_linear)
            _write_fortran_record(f, tau_data)

            # Records 7 to 7+nlevel-1: Departure coefficients (level-major)
            # The input `depart_coeffs` is already (nlevel, ndep), which is correct.
            for j in range(nlevel):
                level_data = depart_coeffs[j, :]
                b_dep_data = struct.pack(f'<{ndep}f', *level_data)
                _write_fortran_record(f, b_dep_data)

    else:
        # Create a formatted ASCII file
        #print(f"Creating ASCII departure file: {filename}")
        with open(filename, 'w') as f:
            # 1. Header (8 lines) - we'll just write placeholders
            f.write("  placeholder_str           0.0         0.0         0.0\n" * 8)

            # 2. Abundance
            f.write(f"{abundance:.2f}\n")

            # 3. Number of depth points
            f.write(f"{ndep}\n")

            # 4. Number of levels
            f.write(f"{nlevel}\n")

            # 5. Optical depths (log10(tau))
            # The format is free, so one value per line is simplest.
            for t in tau_log10:
                f.write(f" {t:12.4f}\n")

            # 6. Departure coefficients (depth-major)
            # We need to transpose the (nlevel, ndep) array to (ndep, nlevel).
            depart_transposed = depart_coeffs.T
            # Use np.savetxt for easy and clean formatting.
            np.savetxt(f, depart_transposed, fmt='%.8f')
#--------------------------------------------------------------------------------

