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
import tempfile
import logging
from astropy.io import ascii
from astropy.table import Table, Column

from ispec.atmospheres import write_atmosphere
from ispec.lines import write_atomic_linelist
from ispec.common import which
from ispec.common import is_moog_support_enabled
from ispec.spectrum import create_spectrum_structure, resample_spectrum
from effects import _filter_linelist, apply_post_fundamental_effects

def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, regions=None, tmp_dir=None, timeout=1800):
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, tmp_dir=tmp_dir, timeout=timeout)

def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, tmp_dir=None, timeout=1800):
    if not is_moog_support_enabled():
        raise Exception("MOOG support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
    moog_dir = ispec_dir + "/synthesizer/moog/"
    moog_executable = moog_dir + "MOOGSILENT"

    #molecules = linelist['molecule'] == "T"
    #if len(np.where(molecules)[0]) > 0:
        #logging.warn("The molecules have been removed from the linelist because if not MOOG is unreasonably slow")
        #linelist = linelist[~molecules]


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

    # MOOG does not have hard-coded lines, we use a external file:
    #hydrogen_lines_file = moog_dir + "/hydrogen_moog_lines.txt"
    #hydrogen_lines = ascii.read(hydrogen_lines_file, names=["wave_A", "spectrum_moog_species", "lower_state_eV", "loggf"])
    hydrogen_lines_file = moog_dir + "DATA/Hlinedata"
    hydrogen_lines = ascii.read(hydrogen_lines_file, names=["wave_A", "spectrum_moog_species", "lower_state_eV", "loggf", "designation"])

    # MOOG is not going to scale the abundances because we are indicating
    # our abundances in the input and that overrides any other prescription, thus
    # we have to manually scale (but do not change Hydrogen and Helium!)
    abundances = abundances.copy()
    efilter = np.logical_and(abundances['code'] != 1, abundances['code'] != 2)
    efilter = np.logical_and(efilter, abundances['code'] <= 92)
    abundances['Abund'][efilter] += MH

    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        for fixed_abundance in fixed_abundances:
            index = np.where(abundances['code'] == fixed_abundance['code'])[0]
            abundances['Abund'][index] = fixed_abundance['Abund']

    synth_fluxes = []
    synth_waveobs = []
    for i, region in enumerate(regions):
        # It is better not to synthesize in a single run a big chunk of wavelength so
        # we split the computation in several pieces
        max_segment = 100. # nm
        if (region['wave_top'] - region['wave_base'])/wave_step > max_segment/wave_step:
            segment_wave_base = np.arange(region['wave_base'], region['wave_top'], max_segment)
            segments = np.recarray((len(segment_wave_base),),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'] = segment_wave_base
            segments['wave_top'] = segment_wave_base + max_segment - wave_step
            segments['wave_top'][-1] = region['wave_top'] # Last segment should not over pass the original global limits
        else:
            segments = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
            segments['wave_base'][0] = region['wave_base']
            segments['wave_top'][0] = region['wave_top']

        for segment in segments:
            wave_base = segment['wave_base']
            wave_top = segment['wave_top']

            tmp_execution_dir = tempfile.mkdtemp(dir=tmp_dir)
            os.makedirs(tmp_execution_dir+"/DATA/")
            atmosphere_filename = tmp_execution_dir + "/model.in"
            linelist_filename = tmp_execution_dir + "/lines.in"
            barklem_linelist_filename = tmp_execution_dir + "/DATA/Barklem.dat"
            stronglinelist_file = tmp_execution_dir + "/stronglines.in"

            if atmosphere_layers_file is None:
                atmosphere_filename = write_atmosphere(atmosphere_layers, teff, logg, MH, code="moog", atmosphere_filename=atmosphere_filename, tmp_dir=tmp_dir)
            else:
                shutil.copyfile(atmosphere_layers_file, atmosphere_filename)

            # Write (MOOG requires a separate file for damping coeff if we want to provide rad coeff and alpha from ABO theory)
            wfilter = np.logical_and(linelist['wave_nm'] >= wave_base, linelist['wave_nm'] <= wave_top)
            linelist_filename = write_atomic_linelist(linelist[wfilter], linelist_filename=linelist_filename, code="moog", tmp_dir=tmp_dir)
            barklem_linelist_filename = write_atomic_linelist(linelist[wfilter], linelist_filename=barklem_linelist_filename, code="moog_barklem", tmp_dir=tmp_dir)

            molecules = linelist['molecule'][wfilter] == 'T'
            num_molecules = len(np.where(molecules)[0])

            # Append microturbulence, solar abundances and metallicity
            moog_atmosphere = open(atmosphere_filename, "a")
            atom_abundances = abundances[np.logical_and(abundances['code'] > 1, abundances['code'] <= 92)] # Don't update hydrogen or helium abundances
            moog_atmosphere.write("  %.2f\n" % (microturbulence_vel))
            moog_atmosphere.write("NATOMS=   %i %.2f\n" % (len(atom_abundances), MH))
            for atom_abundance in atom_abundances:
                abund = 12.036 + atom_abundance['Abund'] # From SPECTRUM format to Turbospectrum
                moog_atmosphere.write("%i  %.2f\n" % (atom_abundance['code'], abund))

            # Molecules are required always, even if it is a linelist without molecules,
            # or MOOG will not compute molecular equilibrium which might affect other internal calculations
            #if num_molecules > 0:
                #unique_molecules = np.unique(linelist['spectrum_moog_species'][wfilter][molecules])
                #moog_atmosphere.write("NMOL      %i\n" % (len(unique_molecules)))
                #for specie in unique_molecules:
                    #moog_atmosphere.write("  %s\n" % (specie))

            # Molecule list as used by Jorge Melendez (private communication)
            moog_atmosphere.write("NMOL      28\n")
            moog_atmosphere.write("  101.0   106.0   107.0   108.0   112.0  126.0\n")
            moog_atmosphere.write("  606.0   607.0   608.0\n")
            moog_atmosphere.write("  707.0   708.0\n")
            moog_atmosphere.write("  808.0   812.0   822.0   823.0   840.0\n")
            moog_atmosphere.write("  10108.0 10820.0 60808.0\n")
            moog_atmosphere.write("  6.1     7.1     8.1   12.1  20.1  22.1  23.1  26.1  40.1\n")
            moog_atmosphere.close()

            # Add hydrogen lines
            # Provide some margin or near-by deep lines might be omitted
            margin = 2. # 2 nm
            wfilter = np.logical_and(hydrogen_lines['wave_A'] >= (wave_base-margin)*10., hydrogen_lines['wave_A'] <= (wave_top+margin)*10.)
            selected_hydrogen_lines = hydrogen_lines[wfilter]
            if len(selected_hydrogen_lines) > 40:
                # TODO: Find a work around to this
                raise Exception("Wavelength range too big, it includes too many hydrogen lines")
            out = open(stronglinelist_file, "w")
            for line in selected_hydrogen_lines:
                out.write("%10.3f%10s%10.3f%10.3f%10s%10s%10s%10s\n" \
                        % (line['wave_A'], line['spectrum_moog_species'], line['lower_state_eV'], line['loggf'], "", "", "", ""))
            out.close()

            par_file = open(tmp_execution_dir + "/batch.par", "w")
            par_file.write("synth\n")
            par_file.write("standard_out moog.std\n")
            par_file.write("summary_out  moog.sum\n")
            par_file.write("smoothed_out moog.spec\n")
            par_file.write("model_in     model.in\n")
            par_file.write("lines_in     lines.in\n")
            par_file.write("stronglines_in stronglines.in\n")
            par_file.write("atmosphere  0\n") # controls the output of atmosphere quantities
            #if num_molecules > 0:
                #par_file.write("molecules   1\n") # controls the molecular equilibrium calculations (1 = do molecular equilibrium but do not print results)
            #else:
                #par_file.write("molecules   0\n")
            # molecules should be always on or moog does not run
            par_file.write("molecules   1\n") # controls the molecular equilibrium calculations (1 = do molecular equilibrium but do not print results)
            par_file.write("lines       1\n") # controls the output of line data (print out standard information about the input line list)
            par_file.write("strong      1\n")
            par_file.write("flux/int    0\n") # choses integrated flux or central intensity (0 = integrated flux calculations)
            par_file.write("damping     1\n") # Use Barklem.dat file with waals, alpha and rad damping coeff, if not found then do like damping = 0
            par_file.write("freeform    0\n") # Linelist format of 7 columns with numbers %10.3f and comment %10s
            par_file.write("plot        3\n")
            par_file.write("abundances  0  1\n")
            par_file.write("isotopes    0  1\n")
            par_file.write("synlimits\n")
            wave_range_of_line_influence = 20.0 # Amstrong
            par_file.write(" %.2f %.2f %.2f %.1f\n" % (wave_base*10., wave_top*10., wave_step*10., wave_range_of_line_influence))
            par_file.write("obspectrum    5\n")
            #par_file.write("plotpars    0\n")
            par_file.write("plotpars    1\n")
            par_file.write(" %.2f %.2f 0.0 1.10\n" % (wave_base*10., wave_top*10.)) # Ploting limits
            par_file.write(" 0.0     0.0    0.00  1.0\n") # vshift       lamshift    obsadd    obsmult
            par_file.write(" r       0.0  0.0  0.0  0.00  0.0") #  smooth-type  FWHM-Gauss  vsini     LimbDarkeningCoeff    FWHM-Macro     FWHM-Loren
            par_file.close()

            previous_cwd = os.getcwd()
            os.chdir(tmp_execution_dir)

            command = moog_executable
            command_input = ""

            # If timeout command exists in PATH, then use it to control moog execution time
            if which("timeout") is not None:
                command = "timeout %i " % (timeout) + command

            # MOOG clears the screen, so it is better to always don't do verbose this part
            # also the information printed is not specially useful
            #if verbose == 1:
                #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
            #else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the moog synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            os.chdir(previous_cwd)

            try:
                data = np.loadtxt(tmp_execution_dir+"/moog.spec", skiprows=2)
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

            shutil.rmtree(tmp_execution_dir)

    # Zero values, when convolved, remain zero so we give a very tiny flux to avoid this problem
    synth_fluxes[synth_fluxes <= 0] = 10e-9

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
