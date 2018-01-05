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
from effects import _filter_linelist, apply_post_fundamental_effects


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, use_molecules=False, tmp_dir=None, timeout=1800):
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)

def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, use_molecules=False, tmp_dir=None, timeout=1800):
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

    return synth_spectrum['flux']


