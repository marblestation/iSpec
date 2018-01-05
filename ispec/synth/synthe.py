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

from ispec.lines import write_atomic_linelist
from ispec.common import which, is_synthe_support_enabled
from ispec.spectrum import create_spectrum_structure, resample_spectrum
from effects import _filter_linelist, apply_post_fundamental_effects


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, molecules_files=None, regions=None, tmp_dir=None, timeout=1800):
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, R=0, macroturbulence=0, vsini=0, limb_darkening_coeff=0, tmp_dir=tmp_dir, timeout=timeout)

def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0,  atmosphere_layers_file=None, linelist_file=None, molecules_files=None, regions=None, R=None, macroturbulence=None, vsini=None, limb_darkening_coeff=None, tmp_dir=None, timeout=1800):
    if not is_synthe_support_enabled():
        raise Exception("Synthe support is not enabled")

    ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
    atmos_dir = ispec_dir + "/synthesizer/atmos/"
    system_64bits = sys.maxsize > 2**32
    if system_64bits:
        xnfpelsyn_executable = atmos_dir + "bin.amd64/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.amd64/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.amd64/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.amd64/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.amd64/synthe.exe"
        spectrv_executable = atmos_dir + "bin.amd64/spectrv.exe"
        rotate_executable = atmos_dir + "bin.amd64/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.amd64/syntoascanga.exe"
    else:
        logging.warn("*************************************************")
        logging.warn("Synthe does not work properly in 32 bits systems!")
        logging.warn("*************************************************")
        xnfpelsyn_executable = atmos_dir + "bin.ia32/xnfpelsyn.exe"
        synbeg_executable = atmos_dir + "bin.ia32/synbeg.exe"
        #rline2.exe # It does not exist in the source code!
        rgfallinesnew_executable = atmos_dir + "bin.ia32/rgfalllinesnew.exe"
        rmolescasc_executable = atmos_dir + "bin.ia32/rmolecasc.exe"
        synthe_executable = atmos_dir + "bin.ia32/synthe.exe"
        spectrv_executable = atmos_dir + "bin.ia32/spectrv.exe"
        rotate_executable = atmos_dir + "bin.ia32/rotate.exe"
        syntoascanga_executable = atmos_dir + "bin.ia32/syntoascanga.exe"
    atmos_molecules = atmos_dir + "lines/molecules.dat"
    atmos_helium = atmos_dir + "lines/he1tables.dat"
    atmos_continua = atmos_dir + "lines/continua.dat"


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


    # Update abundances with the ones that should be fixed to a given value and
    # not affected by metallicity scalation
    if fixed_abundances is not None and len(fixed_abundances) > 0:
        abundances = abundances.copy()
        for fixed_abundance in fixed_abundances:
            index = np.where(abundances['code'] == fixed_abundance['code'])[0]
            abundances['Abund'][index] = fixed_abundance['Abund'] - MH



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
            os.symlink(atmos_molecules, tmp_execution_dir+"/fort.2")
            os.symlink(atmos_helium, tmp_execution_dir+"/fort.18")
            os.symlink(atmos_continua, tmp_execution_dir+"/fort.17")


            previous_cwd = os.getcwd()
            os.chdir(tmp_execution_dir)

            # XNFPELSYN pretabulates continuum opacities and number densities for
            # different chemical elements and writes them to fort.10
            command = xnfpelsyn_executable

            command_input = "SURFACE INTENSI 17 1.,.9,.8,.7,.6,.5,.4,.3,.25,.2,.15,.125,.1,.075,.05,.025,.01\n"
            #command_input = "SURFACE FLUX\n"
            command_input += "ITERATIONS 1 PRINT 2 PUNCH 2\n"
            command_input += "CORRECTION OFF\n"
            command_input += "PRESSURE OFF\n"
            command_input += "READ MOLECULES\n"
            command_input += "MOLECULES ON\n"
            command_input += "TEFF   %.0f  GRAVITY %.5f LTE\n" % (teff, logg)
            command_input += "TITLE  ISPEC\n"
            command_input += " OPACITY IFOP 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 0 0 0 0 0\n"
            mixing_length_param = 1.25
            command_input += " CONVECTION ON   %.2f TURBULENCE OFF  0.00  0.00  0.00  0.00\n" % (mixing_length_param)
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
            model = command_input

            # If timeout command exists in PATH, then use it to control synthe execution time
            which_timeout = which("timeout")
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            # Never verbose because WIDTH's output is printed on stdout (not saved on a file)
            # and it is needed by iSpec
            #if verbose == 1:
                #proc = subprocess.Popen(command.split(), stdin=subprocess.PIPE)
            #else:
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode
            #print out

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            # synbeg: Reads the fundamental parameters of the process: wavelength range,
            # resolution of the final spectrum before degrading, whether unclassified
            # (pre-dicted) lines should be included, etc., and writes them to fort.93
            command = synbeg_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            # microturbulence is added by summing the squares to the microturbulence indicated in the mode atmosphere, which we have forced to zero before
            turbv = microturbulence_vel
            # NOTE: For the wavelength range, it does not work as other codes... it will ignore any lines outside the provided range
            # even if they are strong near-by lines that will affect the region of interest. So we increase the region to synthesize
            # and we will cut later.
            # Provide some margin or near-by deep lines might be omitted
            margin = 2. # 2 nm
            #command_input =  "AIR        %-9.1f %-9.1f 600000. %9.2f    0     30    .0001     1    0\n" % (wave_base-margin, wave_top+margin, turbv)
            command_input =  "AIR        %-9.1f %-9.1f 600000. %9.2f    0     10    .001      0    0\n" % (wave_base-margin, wave_top+margin, turbv)
            command_input += "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD\n"
            #AIR indicates that the wavelengths are in AIR. VAC would provide vacuum wavelengths
            #WLBEG and WLEND are the starting and ending points of the synthesis, in nanometers
            #RESOLU is the resolution at which the calculation is performed. Practically, SYNTHE calculates the transfer through the atmosphere at wavelength intervals with such spacing. Of course, reducing the resolution will lead to a faster calculation, but also to a poorer sampling of the radiative transfer through the atmosphere. We thus suggest not to go below a resolution of 100000. This value is adequate for comparison with high resolution observed spectra.
            #TURBV is the microturbulence we want SYNTHE to add to the one in the atmosphere model. Since microturbulence is added by summing the squares, and we have a VTURB=1 model, we need to add 1.67 to obtain the final 1.95 km/s.
            #IFNLTE is set to 0 because we want a LTE calculation
            #LINOUT if it is negative, line data are not saved and it speeds up the process
            #CUTOFF is used to keep the weakest transitions out of the output files. With this setting, any absorption subtracting at its center less than 1/10000 of the intensity will be cut off.

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode
            #print out

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")


            if linelist_file is None:
                # Provide some margin or near-by deep lines might be omitted
                margin = 2. # 2 nm
                wfilter = np.logical_and(linelist['wave_nm'] >= wave_base-margin, linelist['wave_nm'] <= wave_top+margin)
                linelist_filename, molecules_filenames = write_atomic_linelist(linelist[wfilter], code="synthe", tmp_dir=tmp_execution_dir)
            else:
                linelist_filename = linelist_file
                molecules_filenames = molecules_files

            ## Compute atomic lines
            os.symlink(linelist_filename, tmp_execution_dir+"/fort.11")
            # rgfallinesnew: adds line information and line opacity data by reading
            # the adequate opacity file from fort.11.
            # The line opacity data are sent to fort.12 (and the line identification
            # to fort.14) if the calculation is LTE; it goes to fort.19 (line
            # identification to fort.20) if not.
            command = rgfallinesnew_executable
            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            #print out
            #print "-"*80
            errcode = proc.returncode
            os.remove(tmp_execution_dir+"/fort.11")

            # Compute molecules
            if molecules_filenames is not None:
                for molecules_filename in molecules_filenames:
                    os.symlink(molecules_filename, tmp_execution_dir+"/fort.11")
                    # RMOLEC is the homologue of RGFALLTEST for diatomic molecules
                    # It writes line opacity to fort.12 and line identifications
                    # to fort. 14
                    command = rmolescasc_executable
                    proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
                    # wait for the process to terminate
                    out, err = proc.communicate()
                    errcode = proc.returncode
                    os.remove(tmp_execution_dir+"/fort.11")
                    #print "Filename:", molecules_filename
                    #print out


            # SYNTHE computes line opacity data based on the fundamental parameters of
            # the model (the ones written by SYNBEG to fort.93) and writes them to fort.93
            command = synthe_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            config_data = "0.0       0.        1.        0.        0.        0.        0.        0.\n"
            config_data += "0.\n"
            config_data += "RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW\n"
            config = open(tmp_execution_dir+"/fort.25", "w")
            config.write(config_data)
            config.close()

            # SPECTRV reads from unit 9 the file written by the program SYNTHE and com-
            # putes the continuum opacities and the overall synthetic spectrum
            # (intensities at 17 angles); this spectrum is then written to fort.7
            command = spectrv_executable
            command_input = model

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")



            # ROTATE handles rotational broadening of spectral lines. It also integrates the
            # intensity data given by SPECTRV to compute total flux.
            # The arguments are the number of values of the projected rotational
            # velocity v sin i (first row) and the velocities themselves (second row)
            #
            # Working notes: "It is recommended to always use SURFACE INTENSITY instead of SURFACE FLUX.
            # we can get flux spectra of non-rotating stars with SURFACE INTENSITY and
            # v sin i = 0"
            command = rotate_executable
            intensities_filename = tmp_execution_dir+"/intensities.bin"
            if os.path.exists(tmp_execution_dir+"/fort.7"):
                os.rename(tmp_execution_dir+"/fort.7", intensities_filename)
            else:
                raise Exception("No intensities were calculated, synthesis failed!")
            #os.remove(tmp_execution_dir+"/fort.1")
            os.symlink(intensities_filename, tmp_execution_dir+"/fort.1")
            #command_input = "    1   50\n"
            command_input = "    1\n"
            command_input +="0.\n"

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate(input=command_input)
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            flux_filename = tmp_execution_dir+"/flux.bin"
            os.rename(tmp_execution_dir+"/ROT1", flux_filename)
            os.remove(tmp_execution_dir+"/fort.1")
            os.remove(tmp_execution_dir+"/fort.2")

            #os.remove(tmp_execution_dir+"/fort.3")
            os.symlink(flux_filename, tmp_execution_dir+"/fort.1")
            os.symlink(tmp_execution_dir+"/lines.txt", tmp_execution_dir+"/fort.3")
            os.symlink(tmp_execution_dir+"/spectrum.txt", tmp_execution_dir+"/fort.2")
            os.symlink(tmp_execution_dir+"/dump.txt", tmp_execution_dir+"/fort.4")


            command = syntoascanga_executable

            # If timeout command exists in PATH, then use it to control synthe execution time
            if which_timeout is not None:
                command = "timeout %i " % (timeout) + command

            proc = subprocess.Popen(command.split(), stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            # wait for the process to terminate
            out, err = proc.communicate()
            errcode = proc.returncode

            if errcode == 124: # TIMEOUT
                logging.error("A timeout has occurred in the synthe synthesis process.")
                raise Exception("Timeout: Synthesis failed!")

            try:
                data = np.loadtxt(tmp_execution_dir+"/spectrum.txt")
            except:
                #print out
                sys.stdout.flush()
                raise Exception("Synthesis failed!")
            synth_waveobs_tmp = data[:,0] / 10. # Armstrong to nm
            # NOTE: We provided an artificially bigger wavelength range, so that synthe will consider near-by deep lines
            # Now we correct that and we reduce the wavelength range to the correct one:
            wfilter = np.logical_and(synth_waveobs_tmp >= wave_base, synth_waveobs_tmp <= wave_top)
            synth_waveobs = np.hstack((synth_waveobs, synth_waveobs_tmp[wfilter]))

            synth_fluxes_tmp = data[:,3]
            synth_fluxes = np.hstack((synth_fluxes, synth_fluxes_tmp[wfilter]))

            os.chdir(previous_cwd)
            shutil.rmtree(tmp_execution_dir)

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




