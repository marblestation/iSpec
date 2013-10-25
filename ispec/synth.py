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
import time
from datetime import datetime, timedelta
import numpy as np
from mpfitmodels import MPFitModel
#from continuum import fit_continuum
from abundances import write_SPECTRUM_abundances, write_SPECTRUM_fixed_abundances, determine_abundances
from atmospheres import write_atmosphere, interpolate_atmosphere_layers
from lines import write_SPECTRUM_linelist
from spectrum import create_spectrum_structure, convolve_spectrum, correct_velocity, resample_spectrum
from multiprocessing import Process
from multiprocessing import Queue
from multiprocessing import JoinableQueue
from Queue import Empty

import tempfile
import log
import logging


def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=900, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, regions=None, waveobs_mask=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs.
    In case regions is specified (recarray with 'wave_base' and 'wave_top'),
    only those regions will be computed.

    No macroturbulence, rotation (vsini), limb darkening coefficient or resolution is considered
    in this process. That's why it is named as "fundamental" spectrum.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be set to 'None'.
    """
    return generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence = 0.0, vsini = 0.0, limb_darkening_coeff = 0.0, R=0, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, regions=regions)


def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, gui_queue=None, timeout=900, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, regions=None, waveobs_mask=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs.
    In case regions is specified (recarray with 'wave_base' and 'wave_top'),
    only those regions will be computed.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be set to 'None'.
    """
    if fixed_abundances is None:
        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])

    # All synthetic spectra contain a -0.28 shift (unknown reason) so we correct it
    # by moving the waveobs in the other sense:
    # ** The shift has been validated comparing with radial velocities from HARPS
    #    for a wide range of samples in the parameter space (benchmark stars)
    rv_shift = 0.28
    spectrum = create_spectrum_structure(waveobs)
    spectrum = correct_velocity(spectrum, rv_shift)
    waveobs = spectrum['waveobs']

    if waveobs_mask is None:
        if regions is None:
            waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
            # Limit linelist
            wave_base = np.min(waveobs)
            wave_top = np.max(waveobs)
            lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
            linelist = linelist[lfilter]
        else:
            waveobs_mask = __create_waveobs_mask(waveobs, regions)
            # Limit linelist
            linelist = __filter_linelist(linelist, regions)

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_fixed_abund_file = False
    remove_tmp_linelist_file = False
    if atmosphere_layers_file is None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
        remove_tmp_atm_file = True
    if abundances_file is None:
        abundances_file = write_SPECTRUM_abundances(abundances)
        remove_tmp_abund_file = True
    if fixed_abundances_file is None:
        fixed_abundances_file = write_SPECTRUM_fixed_abundances(fixed_abundances)
        remove_tmp_fixed_abund_file = True
    if linelist_file is None:
        linelist_file = write_SPECTRUM_linelist(linelist)
        remove_tmp_linelist_file = True
    nlayers = len(atmosphere_layers)

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    #process_communication_queue = Queue()
    process_communication_queue = JoinableQueue()

    p = Process(target=__generate_spectrum, args=(process_communication_queue, waveobs, waveobs_mask, atmosphere_layers_file, linelist_file, abundances_file, fixed_abundances_file), kwargs={'microturbulence_vel': microturbulence_vel, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose})
    p.start()
    fluxes = np.zeros(len(waveobs))
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) == np.ndarray:
                # Results received!
                fluxes = data
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
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the synthetic spectrum generation.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("The synthetic spectrum generation has failed for these astrophysical parameters.")
        p.terminate()
    else:
        p.join()

    if remove_tmp_atm_file:
        os.remove(atmosphere_layers_file)
    if remove_tmp_abund_file:
        os.remove(abundances_file)
    if remove_tmp_fixed_abund_file:
        os.remove(fixed_abundances_file)
    if remove_tmp_linelist_file:
        os.remove(linelist_file)

    return fluxes

def __enqueue_progress(process_communication_queue, v):
    process_communication_queue.put(("self.update_progress(%i)" % v))
    process_communication_queue.join()

def __generate_spectrum(process_communication_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, nlayers=56, verbose=0):
    """
    Generate synthetic spectrum and apply macroturbulence, rotation (visini), limb darkening coeff and resolution except
    if all those parameters are set to zero, in that case the fundamental synthetic spectrum is returned.
    """
    import synthesizer

    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, 0, nlayers, verbose, update_progress_func)
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10
    if R > 0:
        # Use iSpec convolution routine instead of SPECTRUM one, since iSpec is more reliable
        spectrum = create_spectrum_structure(waveobs, fluxes)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']
    process_communication_queue.put(fluxes)



def apply_post_fundamental_effects(waveobs, fluxes, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, gui_queue=None, timeout=900):
    """
    Apply macroturbulence, rotation (vsini), limb darkening coefficient and/or resolution
    """

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    process_communication_queue = Queue()

    p = Process(target=__apply_post_fundamental_effects, args=(process_communication_queue, waveobs, fluxes), kwargs={'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'verbose': verbose})
    p.start()
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            data = process_communication_queue.get(timeout=1)
            if type(data) == np.ndarray:
                # Results received!
                fluxes = data
                break
            elif gui_queue is not None:
                # GUI update
                # It allows communications between process in order to update the GUI progress bar
                gui_queue.put(data)
        except Empty:
            # No results, continue waiting
            pass
        num_seconds += 1
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the application of post fundamental effects.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("The application of post fundamental effects has failed.")
        p.terminate()
    else:
        p.join()

    return fluxes


def __apply_post_fundamental_effects(process_communication_queue, waveobs, fluxes, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0):
    """
    Apply macroturbulence, rotation (visini), limb darkening coeff and resolution to already generated fundamental synthetic spectrum.
    """
    import synthesizer
    #update_progress_func = lambda v: process_communication_queue.put(("self.update_progress(%i)" % v))
    update_progress_func = lambda v: __enqueue_progress(process_communication_queue, v)
    fluxes = synthesizer.apply_post_fundamental_effects(waveobs*10., fluxes, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, 0, verbose, update_progress_func)
    # Avoid zero fluxes, set a minimum value so that when it is convolved it
    # changes. This way we reduce the impact of the following problem:
    # SPECTRUM + MARCS makes some strong lines to have zero fluxes (i.e. 854.21nm)
    zeros = np.where(fluxes <= 1.0e-10)[0]
    fluxes[zeros] = 1.0e-10
    if R > 0:
        # Use iSpec convolution routine instead of SPECTRUM one, since iSpec is more reliable
        spectrum = create_spectrum_structure(waveobs, fluxes)
        fluxes = convolve_spectrum(spectrum, R, from_resolution=None, frame=None)['flux']
    process_communication_queue.put(fluxes)


class SynthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, abundances, teff=5000, logg=3.0, MH=0.0, vmic=2.0, vmac=0.0, vsini=2.0, limb_darkening_coeff=0.0, R=0, continuum_correction=1.0):
        self.elements = {}
        #self.elements["1"] = "H"
        #self.elements["2"] = "He"
        self.elements["3"] = "Li"
        self.elements["4"] = "Be"
        self.elements["5"] = "B"
        self.elements["6"] = "C"
        self.elements["7"] = "N"
        self.elements["8"] = "O"
        self.elements["9"] = "F"
        self.elements["10"] = "Ne"
        self.elements["11"] = "Na"
        self.elements["12"] = "Mg"
        self.elements["13"] = "Al"
        self.elements["14"] = "Si"
        self.elements["15"] = "P"
        self.elements["16"] = "S"
        self.elements["17"] = "Cl"
        self.elements["18"] = "Ar"
        self.elements["19"] = "K"
        self.elements["20"] = "Ca"
        self.elements["21"] = "Sc"
        self.elements["22"] = "Ti"
        self.elements["23"] = "V"
        self.elements["24"] = "Cr"
        self.elements["25"] = "Mn"
        self.elements["26"] = "Fe"
        self.elements["27"] = "Co"
        self.elements["28"] = "Ni"
        self.elements["29"] = "Cu"
        self.elements["30"] = "Zn"
        self.elements["31"] = "Ga"
        self.elements["32"] = "Ge"
        self.elements["33"] = "As"
        self.elements["34"] = "Se"
        self.elements["35"] = "Br"
        self.elements["36"] = "Kr"
        self.elements["37"] = "Rb"
        self.elements["38"] = "Sr"
        self.elements["39"] = "Y"
        self.elements["40"] = "Zr"
        self.elements["41"] = "Nb"
        self.elements["42"] = "Mo"
        self.elements["43"] = "Tc"
        self.elements["44"] = "Ru"
        self.elements["45"] = "Rh"
        self.elements["46"] = "Pd"
        self.elements["47"] = "Ag"
        self.elements["48"] = "Cd"
        self.elements["49"] = "In"
        self.elements["50"] = "Sn"
        self.elements["51"] = "Sb"
        self.elements["52"] = "Te"
        self.elements["53"] = "I"
        self.elements["54"] = "Xe"
        self.elements["55"] = "Cs"
        self.elements["56"] = "Ba"
        self.elements["57"] = "La"
        self.elements["58"] = "Ce"
        self.elements["59"] = "Pr"
        self.elements["60"] = "Nd"
        self.elements["61"] = "Pm"
        self.elements["62"] = "Sm"
        self.elements["63"] = "Eu"
        self.elements["64"] = "Gd"
        self.elements["65"] = "Tb"
        self.elements["66"] = "Dy"
        self.elements["67"] = "Ho"
        self.elements["68"] = "Er"
        self.elements["69"] = "Tm"
        self.elements["70"] = "Yb"
        self.elements["71"] = "Lu"
        self.elements["72"] = "Hf"
        self.elements["73"] = "Ta"
        self.elements["74"] = "W"
        self.elements["75"] = "Re"
        self.elements["76"] = "Os"
        self.elements["77"] = "Ir"
        self.elements["78"] = "Pt"
        self.elements["79"] = "Au"
        self.elements["80"] = "Hg"
        self.elements["81"] = "Tl"
        self.elements["82"] = "Pb"
        self.elements["83"] = "Bi"
        self.elements["84"] = "Po"
        self.elements["85"] = "At"
        self.elements["86"] = "Rn"
        self.elements["87"] = "Fr"
        self.elements["88"] = "Ra"
        self.elements["89"] = "Ac"
        self.elements["90"] = "Th"
        self.elements["91"] = "Pa"
        self.elements["92"] = "U"
        self.elements["101"] = "Md"
        self.elements["106"] = "Sg"
        self.elements["107"] = "Bh"
        self.elements["108"] = "Hs"
        self.elements["112"] = "Cn"
        self.elements["113"] = "Uut"
        self.elements["114"] = "Uuq"

        self.modeled_layers_pack = modeled_layers_pack
        self.linelist = linelist
        self.abundances = abundances
        #
        self.calculation_time = 0
        self.waveobs = None
        self.waveobs_mask = None
        self.cache = {}
        p = [teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, R, continuum_correction]
        #
        self.abundances_file = None
        self.linelist_file = None
        self.noise = None
        self.last_continuum_correction = None
        super(SynthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        # Consider new abundances as fixed
        fixed_abundances = self.free_abundances()

        abundances_key = " ".join(map(str, fixed_abundances['Abund']))
        complete_key = "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %i %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()), self.continuum_correction())
        complete_key += abundances_key
        key = "%.2f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
        key += abundances_key
        if self.cache.has_key(key):
            print "Cache:", complete_key
            self.last_fluxes = self.cache[key].copy()
        else:
            print "Generating:", complete_key


            # Atmosphere
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH())
            # Fundamental synthetic fluxes
            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, atmosphere_layers, self.teff(), self.logg(), self.MH(), self.linelist, self.abundances, fixed_abundances, microturbulence_vel=self.vmic(), abundances_file=self.abundances_file, linelist_file=self.linelist_file, waveobs_mask=self.waveobs_mask, verbose=0)

            # Fit continuum
            #if self.fit_continuum_func is not None:
                #spectrum = create_spectrum_structure(self.waveobs)
                #spectrum['flux'] = self.last_fluxes
                ## Add noise
                #if self.noise is not None:
                    #spectrum['flux'] += self.noise
                #continuum_model = self.fit_continuum_func(spectrum)
                #self.last_continuum_correction = continuum_model(self.waveobs)
                #inormalize = np.where(self.last_continuum_correction != 0)[0]
                #self.last_fluxes[inormalize] /= self.last_continuum_correction[inormalize]
            # Optimization to avoid too small changes in parameters or repetition
            self.cache[key] = self.last_fluxes.copy()

        self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=self.R(), verbose=0)
        self.last_final_fluxes *= self.continuum_correction()


        #self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_final_fluxes, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=self.R(), verbose=0)

        return self.last_final_fluxes[self.comparing_mask]

    def fitData(self, waveobs, waveobs_mask, comparing_mask, fluxes, weights=None, parinfo=None, max_iterations=20, quiet=True, fit_continuum_func=None, noise=None):
        self.noise = noise
        self.fit_continuum_func = fit_continuum_func
        base = 9
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.waveobs = waveobs
        self.waveobs_mask = waveobs_mask # Synthesis for wavelengths with mask different from 0.0
        self.comparing_mask = comparing_mask == 1.0 # Wavelengths to be compared for the least square algorithm
        if weights is None:
            weights = np.ones(len(waveobs))
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 1.0   # Residuals are limited between -1.0 and 1.0 (np.tanh(residuals/1.0)) minimizing the influence of bad fluxes
                     # * Spectrum must be a normalized one
        _t0 = default_timer()

        # Write abundances and linelist to avoid writing the same info in each iteration
        self.abundances_file = write_SPECTRUM_abundances(self.abundances)
        self.linelist_file = write_SPECTRUM_linelist(self.linelist)

        super(SynthModel, self).fitData(waveobs[self.comparing_mask], fluxes[self.comparing_mask], weights=weights[self.comparing_mask], parinfo=parinfo, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet)

        residuals = self.last_final_fluxes[self.comparing_mask] - fluxes[self.comparing_mask]
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
        # Chisq using tanh
        self.chisq = np.sum(np.tanh(weights[self.comparing_mask] * residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof
        # Chisq without using tanh for minimizing outliers
        self.basic_chisq = np.sum((weights[self.comparing_mask] * residuals)**2)
        self.reduced_basic_chisq = self.basic_chisq / self.m.dof

        self.cache = {}

        os.remove(self.abundances_file)
        os.remove(self.linelist_file)
        self.abundances_file = None
        self.linelist_file = None

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self): return self._parinfo[0]['value']
    def logg(self): return self._parinfo[1]['value']
    def MH(self): return self._parinfo[2]['value']
    def vmic(self): return self._parinfo[3]['value']
    def vmac(self): return self._parinfo[4]['value']
    def vsini(self): return self._parinfo[5]['value']
    def limb_darkening_coeff(self): return self._parinfo[6]['value']
    def R(self): return self._parinfo[7]['value']
    def continuum_correction(self): return self._parinfo[8]['value']
    def free_abundances(self):
        base = 9
        fixed_abundances = np.recarray((len(self._parinfo)-base, ), dtype=[('code', int),('Abund', float)])
        for i in xrange(len(self._parinfo)-base):
            fixed_abundances['code'] = int(self._parinfo[base+i]['parname'])
            fixed_abundances['Abund'] = self._parinfo[base+i]['value']
        return fixed_abundances

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def eMH(self): return self.m.perror[2]
    def evmic(self): return self.m.perror[3]
    def evmac(self): return self.m.perror[4]
    def evsini(self): return self.m.perror[5]
    def elimb_darkening_coeff(self): return self.m.perror[6]
    def eR(self): return self.m.perror[7]
    def econtinuum_correction(self): return self.m.perror[8]
    def efree_abundances(self):
        base = 9
        eabundances = []
        for i in xrange(len(self._parinfo)-base):
            eabundances.append(self.m.perror[base+i])
        return eabundances

    def transformed_free_abundances(self):
        free_abundances = self.free_abundances()
        efree_abundances = self.efree_abundances()

        transformed_abund = np.recarray((len(free_abundances), ), dtype=[('code', int),('Abund', float), ('element', '|S5'), ('[X/H]', float), ('A(X)', float), ('[X/Fe]', float), ('eAbund', float), ('e[X/H]', float), ('e[X/Fe]', float), ('eA(X)', float)])

        #function abundances, sme, solar_relative
        #abund = sme.abund
        #solar_abund, solar, labels

        #solar_relative = abund-solar

        #abund[2:*] += sme.feh ; rescale by metallicity
        #abund[1:*] = 10^abund[1:*] ; transform to linear
        #return, alog10(abund / abund[0]) + 12
        #end
        sun_log_Nh_over_Ntotal = self.abundances['Abund'][self.abundances['code'] == 1]
        for i in xrange(len(free_abundances)):
            sun_log_Nx_over_Ntotal = self.abundances['Abund'][self.abundances['code'] == free_abundances['code'][i]]
            x_absolute = free_abundances['Abund'][i] + 12. - sun_log_Nh_over_Ntotal # absolute, A(X)
            #x_over_fe = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
            #x_over_h = x_over_fe + self.MH()
            x_over_h = free_abundances['Abund'][i] - sun_log_Nx_over_Ntotal
            x_over_fe = x_over_h - self.MH()

            element = self.elements[str(free_abundances['code'][i])]

            transformed_abund['code'][i] = free_abundances['code'][i]
            transformed_abund['Abund'][i] = free_abundances['Abund'][i]
            transformed_abund['element'][i] = element
            transformed_abund['[X/H]'][i] = x_over_h
            transformed_abund['[X/Fe]'][i] = x_over_fe
            transformed_abund['A(X)'][i] = x_absolute
            transformed_abund['eAbund'][i] = efree_abundances[i]
            transformed_abund['e[X/H]'][i] = efree_abundances[i]
            transformed_abund['e[X/Fe]'][i] = efree_abundances[i]
            transformed_abund['eA(X)'][i] = efree_abundances[i]
        return transformed_abund

    def print_solution(self):
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","vmic","vmac","vsini","limb","R","Cont")
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i\t%8.3f" % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()), self.continuum_correction())
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i\t%8.3f" % (self.eteff(), self.elogg(), self.eMH(), self.evmic(), self.evmac(), self.evsini(), self.elimb_darkening_coeff(), int(self.eR()), self.econtinuum_correction())

        # Append free individual abundances
        abundances_header = ""
        abundances_solution = ""
        abundances_errors = ""

        transformed_abund = self.transformed_free_abundances()
        for i in xrange(len(transformed_abund)):
            element = transformed_abund['element'][i]
            x_absolute_name = "A(" + element + ")"
            x_over_h_name = "[" + element + "/H]"
            x_over_fe_name = "[" + element + "/Fe]"
            x = transformed_abund['Abund'][i]
            x_absolute = transformed_abund['A(X)'][i]
            x_over_h = transformed_abund['[X/H]'][i]
            x_over_fe = transformed_abund['[X/Fe]'][i]
            ex = transformed_abund['eAbund'][i]
            ex_absolute = transformed_abund['eA(X)'][i]
            ex_over_h = transformed_abund['e[X/H]'][i]
            ex_over_fe = transformed_abund['e[X/Fe]'][i]
            abundances_header += "%8s\t%8s\t%8s\t%8s" % (element, x_absolute_name, x_over_h_name, x_over_fe_name)
            abundances_solution += "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (x, x_absolute, x_over_h, x_over_fe)
            abundances_errors += "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (ex, ex_absolute, ex_over_h, ex_over_fe)

        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""
        if len(transformed_abund) > 0:
            print "           ", abundances_header
            print "Abundances:", abundances_solution
            print "Ab. errors:", abundances_errors
            print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","chisq(tanh)","rchisq(tanh)","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.chisq, self.reduced_chisq, self.basic_chisq, self.reduced_basic_chisq, self.rms)
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status

def __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_continuum_correction, free_params, free_abundances, teff_range, logg_range, MH_range):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 9
    free_params = [param.lower() for param in free_params]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base+len(free_abundances))]
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = 100.0 # For auto-derivatives
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [np.min(teff_range), np.max(teff_range)]
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = 0.10 # For auto-derivatives
    #parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [np.min(logg_range), np.max(logg_range)]
    #
    parinfo[2]['parname'] = "MH"
    parinfo[2]['value'] = initial_MH
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = 0.05 # For auto-derivatives
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [np.min(MH_range), np.max(MH_range)]
    #
    parinfo[3]['parname'] = "Vmic"
    parinfo[3]['value'] = initial_vmic
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = 0.5 # For auto-derivatives
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [0.0, 50.0]
    #
    parinfo[4]['parname'] = "Vmac"
    parinfo[4]['value'] = initial_vmac
    parinfo[4]['fixed'] = not parinfo[4]['parname'].lower() in free_params
    parinfo[4]['step'] = 2.00 # For auto-derivatives
    parinfo[4]['limited'] = [True, True]
    parinfo[4]['limits'] = [0.0, 50.0]
    #
    parinfo[5]['parname'] = "Vsini"
    parinfo[5]['value'] = initial_vsini
    parinfo[5]['fixed'] = not parinfo[5]['parname'].lower() in free_params
    parinfo[5]['step'] = 0.1 # For auto-derivatives
    parinfo[5]['limited'] = [True, True]
    parinfo[5]['limits'] = [0.0, 50.0]
    #
    parinfo[6]['parname'] = "limb_darkening_coeff"
    parinfo[6]['value'] = initial_limb_darkening_coeff
    parinfo[6]['fixed'] = not parinfo[6]['parname'].lower() in free_params
    parinfo[6]['step'] = 0.20 # For auto-derivatives
    parinfo[6]['limited'] = [True, True]
    parinfo[6]['limits'] = [0.0, 1.0]
    #
    parinfo[7]['parname'] = "R"
    parinfo[7]['value'] = initial_R
    parinfo[7]['fixed'] = not parinfo[7]['parname'].lower() in free_params
    parinfo[7]['step'] = 100.0 # For auto-derivatives
    parinfo[7]['limited'] = [True, True]
    parinfo[7]['limits'] = [500.0, 300000.0]
    #
    parinfo[8]['parname'] = "continuum_correction"
    parinfo[8]['value'] = initial_continuum_correction
    parinfo[8]['fixed'] = not parinfo[8]['parname'].lower() in free_params
    parinfo[8]['step'] = 0.001 # For auto-derivatives
    parinfo[8]['limited'] = [True, True]
    parinfo[8]['limits'] = [0.98, 1.02]
    #
    base = 9
    for i in xrange(len(free_abundances)):
        parinfo[base+i]['parname'] = str(free_abundances['code'][i])
        parinfo[base+i]['value'] = free_abundances['Abund'][i]
        parinfo[base+i]['fixed'] = not parinfo[base+i]['parname'].lower() in free_params
        parinfo[base+i]['step'] = 0.05 # For auto-derivatives
        parinfo[base+i]['limited'] = [True, True]
        parinfo[base+i]['limits'] = [-30., 0.]

    return parinfo

def __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range):
    """
    Creates the structure needed for the mpfitmodel
    """
    base = 4
    #free_params = ["teff", "logg", "vmic", "mh"]
    free_params = ["teff", "logg", "vmic"]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(base)]
    #
    parinfo[0]['parname'] = "teff"
    parinfo[0]['value'] = initial_teff
    parinfo[0]['fixed'] = not parinfo[0]['parname'].lower() in free_params
    parinfo[0]['step'] = 500.0 # For auto-derivatives
    #parinfo[0]['mpside'] = 2
    #parinfo[0]['mpmaxstep'] = parinfo[0]['step'] * 1.5
    parinfo[0]['limited'] = [True, True]
    parinfo[0]['limits'] = [np.min(teff_range), np.max(teff_range)]
    #
    parinfo[1]['parname'] = "logg"
    parinfo[1]['value'] = initial_logg
    parinfo[1]['fixed'] = not parinfo[1]['parname'].lower() in free_params
    parinfo[1]['step'] = 0.5 # For auto-derivatives
    #parinfo[1]['mpside'] = 2
    parinfo[1]['mpmaxstep'] = 0.50 # Maximum change to be made in the parameter
    #parinfo[1]['mpmaxstep'] = parinfo[1]['step'] * 1.5
    parinfo[1]['limited'] = [True, True]
    parinfo[1]['limits'] = [np.min(logg_range), np.max(logg_range)]
    #
    parinfo[2]['parname'] = "Vmic"
    parinfo[2]['value'] = initial_vmic
    parinfo[2]['fixed'] = not parinfo[2]['parname'].lower() in free_params
    parinfo[2]['step'] = 0.50 # For auto-derivatives
    #parinfo[2]['mpside'] = 2
    #parinfo[2]['mpmaxstep'] = parinfo[2]['step'] * 2.0
    parinfo[2]['limited'] = [True, True]
    parinfo[2]['limits'] = [0., 50.0]
    #
    parinfo[3]['parname'] = "MH"
    parinfo[3]['value'] = initial_MH
    parinfo[3]['fixed'] = not parinfo[3]['parname'].lower() in free_params
    parinfo[3]['step'] = 0.05 # For auto-derivatives
    #parinfo[3]['mpside'] = 2
    #if not parinfo[3]['fixed']:
        #parinfo[3]['mpmaxstep'] = parinfo[3]['step'] * 1.5
    parinfo[3]['limited'] = [True, True]
    parinfo[3]['limits'] = [np.min(MH_range), np.max(MH_range)]

    return parinfo


def __filter_linelist(linelist, segments):
    # Build wavelength points from regions
    lfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if lfilter is None:
            lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
        else:
            lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.))

    if lfilter is not None:
        return linelist[lfilter]
    else:
        return linelist

def __create_waveobs_mask(waveobs, segments):
    # Build wavelength points from regions
    wfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if wfilter is None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_mask = np.zeros(len(waveobs))
    waveobs_mask[wfilter] = 1.0 # Compute fluxes only for selected segments

    return waveobs_mask

def __filter_linemasks_not_in_segments(linemasks, segments):
    if segments is None:
        return linemasks
    else:
        lfilter = linemasks['wave_base'] == -1
        for i, region in enumerate(linemasks):
            wave_base = region['wave_base']
            wave_top = region['wave_top']

            # Consider only lines that are inside segments
            in_segment1 = np.logical_and(segments['wave_base'] <= wave_base, segments['wave_top'] >= wave_base)
            in_segment2 = np.logical_and(segments['wave_base'] <= wave_top, segments['wave_top'] >= wave_top)
            in_segment = np.logical_and(in_segment1, in_segment2)
            if np.all(in_segment == False):
                lfilter[i] = False
            else:
                lfilter[i] = True

        return linemasks[lfilter]

def __create_comparing_mask(waveobs, linemasks, segments):
    # Build wavelength points from regions
    wfilter = None
    for region in linemasks:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        # Consider only lines that are inside segments
        if segments is not None:
            in_segment1 = np.logical_and(segments['wave_base'] <= wave_base, segments['wave_top'] >= wave_base)
            in_segment2 = np.logical_and(segments['wave_base'] <= wave_top, segments['wave_top'] >= wave_top)
            in_segment = np.logical_and(in_segment1, in_segment2)
            if np.all(in_segment == False):
                continue

        if wfilter is None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_linemask = np.zeros(len(waveobs))
    waveobs_linemask[wfilter] = 1.0 # Consider fluxes only for selected line masks

    return waveobs_linemask

def __get_stats_per_linemask(waveobs, fluxes, synthetic_fluxes, weights, free_params, linemasks, verbose=False):

    results = np.recarray((len(linemasks), ), dtype=[('wave_peak', float),('wave_base', float),('wave_top', float),('chisq(tanh)', float),('rchisq(tanh)', float),('chisq', float),('rchisq', float),('rms', float)])
    results['wave_peak'] = linemasks['wave_peak']
    results['wave_base'] = linemasks['wave_base']
    results['wave_top'] = linemasks['wave_top']

    i = 0
    for region in linemasks:
        wave_peak = region['wave_peak']
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        # Do not compare negative or zero fluxes
        wfilter = np.logical_and(wfilter, fluxes > 0.0)

        # Degrees of freedom
        dof = len(waveobs[wfilter]) - len(free_params)
        if dof > 0:
            residuals = synthetic_fluxes[wfilter] - fluxes[wfilter]
            rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
            # Chisq without using tanh
            chisq = np.sum(np.tanh(weights[wfilter] * residuals)**2)
            reduced_chisq = chisq / dof
            # Chisq without using tanh for minimizing outliers
            basic_chisq = np.sum((weights[wfilter] * residuals)**2)
            reduced_basic_chisq = basic_chisq / dof
        else:
            rms = -9999
            chisq = -9999
            reduced_chisq = -9999
            basic_chisq = -9999
            reduced_basic_chisq = -9999

        results['rms'][i] = rms
        results['chisq(tanh)'][i] = chisq
        results['rchisq(tanh)'][i] = reduced_chisq
        results['chisq'][i] = basic_chisq
        results['rchisq'][i] = reduced_basic_chisq
        if verbose:
            header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("wave_peak","wave_base","wave_top","chisq(tanh)","rchisq(tanh)","chisq","rchisq","rms")
            stats = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (wave_peak, wave_base, wave_top, chisq, reduced_chisq, basic_chisq, reduced_basic_chisq, rms)
            if i == 0:
                print "         ", header
            print "Line     ", stats
        i += 1

    return results

def modelize_spectrum(spectrum, continuum_model, modeled_layers_pack, linelist, abundances, free_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, segments=None, linemasks=None, max_iterations=20):
    """
    It matches synthetic spectrum to observed spectrum by applying a least
    square algorithm.

    - free_params is an array that can contain any combination of the following
      strings: ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    - free_abundances can be set to 'None'
    - If segments are specified, the synthetic spectrum will be only generated for
      those regions.
    - If linemasks are specified, only those regions will be used for comparison.
    - It does not compare negative or zero fluxes
    """
    # Normalize
    spectrum_orig = spectrum
    if segments is not None:
        # Build wavelength points from regions
        wfilter = None
        for region in segments:
            wave_base = region["wave_base"]
            wave_top = region["wave_top"]

            if wfilter is None:
                wfilter = np.logical_and(spectrum_orig['waveobs'] >= wave_base, spectrum_orig['waveobs'] <= wave_top)
            else:
                wfilter = np.logical_or(wfilter, np.logical_and(spectrum_orig['waveobs'] >= wave_base, spectrum_orig['waveobs'] <= wave_top))
        spectrum = create_spectrum_structure(spectrum_orig['waveobs'][wfilter], spectrum_orig['flux'][wfilter], spectrum_orig['err'][wfilter])
    else:
        spectrum = create_spectrum_structure(spectrum_orig['waveobs'], spectrum_orig['flux'], spectrum_orig['err'])

    # Normalization
    continuum = continuum_model(spectrum['waveobs'])
    inormalize = np.where(continuum != 0)[0]
    spectrum['flux'][inormalize] /= continuum[inormalize]
    spectrum['err'][inormalize] /= continuum[inormalize]

    waveobs = spectrum['waveobs']
    flux = spectrum['flux']
    err = spectrum['err']


    if free_abundances is None:
        # No fixed abundances
        free_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])

    if segments is None:
        waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
        # Limit linelist
        wave_base = np.min(waveobs)
        wave_top = np.max(waveobs)
        lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
        linelist = linelist[lfilter]
    else:
        waveobs_mask = __create_waveobs_mask(waveobs, segments)
        # Limit linelist
        linelist = __filter_linelist(linelist, segments)

    if linemasks is None:
        comparing_mask = np.ones(len(waveobs)) # Compare all fluxes
    else:
        linemasks = __filter_linemasks_not_in_segments(linemasks, segments)
        # Compare fluxes inside line masks that belong to a segment
        comparing_mask = __create_comparing_mask(waveobs, linemasks, segments)
    # Do not compare negative or zero fluxes
    negative_zero = flux <= 0.0
    comparing_mask[negative_zero] = 0.0

    weights = np.ones(len(waveobs))

    ### Use errors:
    ## Do not compare negative or zero errors
    #negative_zero = err <= 0.0
    #comparing_mask[negative_zero] = 0.0
    #weights[comparing_mask == 1.0] = 1. / spectrum['err'][comparing_mask == 1.0]

    ### Estimate from SNR
    #snr = 100
    #weights[comparing_mask == 1.0] = 1. / (flux[comparing_mask == 1.0] / snr)

    ### Prioritize peaks
    # In the following situation, we want to prioritize synthetic spectra that better reproduce the lines' peaks
    # than the upper part (more uncertainties due to continuum normalization or blending). So we would like to
    # have a lower chisq in the case 35 than in case 34:
    #
    #In [32]: flux = np.asarray([1., 0.8, 0.5, 0.8, 1.0])]
    #
    #In [33]: np.sum(np.tanh((1. / (flux * 0.10)) * (flux - np.asarray([0.9, 0.7, 0.5, 0.7, 0.9])))**2)
    #Out[33]: 2.5992215844110831
    #
    #In [34]: np.sum(np.tanh((1. / (flux * 0.10)) * (flux - np.asarray([0.9, 0.7, 0.51, 0.7, 0.9])))**2)
    #Out[34]: 2.6381786014449662
    #
    #In [35]: np.sum(np.tanh((1. / (flux * 0.10)) * (flux - np.asarray([0.91, 0.7, 0.5, 0.7, 0.9])))**2)
    #Out[35]: 2.5322785648767674
    #weights[comparing_mask == 1.0] = 1. / (flux[comparing_mask == 1.0] * 0.10)
    #weights[comparing_mask == 1.0] = 1. / (flux[comparing_mask == 1.0] * 0.01)
    #weights[comparing_mask == 1.0] = 1. / (flux[comparing_mask == 1.0] * 0.0025)


    teff_range = modeled_layers_pack[3]
    logg_range = modeled_layers_pack[4]
    MH_range = modeled_layers_pack[5]

    initial_continuum_correction = 1.0 # TODO: Remove everywhere the continuum_correction parameters because it is not useful
    parinfo = __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_continuum_correction, free_params, free_abundances, teff_range, logg_range, MH_range)

    synth_model = SynthModel(modeled_layers_pack, linelist, abundances)

    #wave_step = 0.001
    #xaxis = np.arange(np.min(spectrum_orig["waveobs"]), np.max(spectrum_orig["waveobs"]), wave_step)
    #resampled_spectrum = resample_spectrum(spectrum_orig, xaxis)
    #snr = estimate_snr(resampled_spectrum['flux'], num_points=10)
    ##snr = 100
    #sigma = flux/snr
    ## Control sigma zero or negative due to zero or negative fluxes
    #noise = np.zeros(len(flux))
    #fnoise = np.where(sigma > 0.0)[0]
    #noise[fnoise] = np.random.normal(0, sigma[fnoise], len(flux[fnoise]))
    noise = None
    #nknots = 1
    ##median_wave_range = 0.05
    #median_wave_range = 0.
    #max_wave_range = 0.1
    #fixed_value = None
    #model = "Polynomy"
    ##fixed_value = 1.0
    ##model = 'Fixed value'
    #fit_continuum_func = lambda x: fit_continuum(x, independent_regions=segments, nknots=nknots, median_wave_range=median_wave_range, max_wave_range=max_wave_range, fixed_value=fixed_value, model=model)
    fit_continuum_func = None
    synth_model.fitData(waveobs, waveobs_mask, comparing_mask, flux, weights=weights, parinfo=parinfo, max_iterations=max_iterations, quiet=False, fit_continuum_func=fit_continuum_func, noise=noise)
    print "\n"
    stats_linemasks = __get_stats_per_linemask(waveobs, flux, synth_model.last_final_fluxes, weights, free_params, linemasks, verbose=True)
    print "\n"
    synth_model.print_solution()

    # Collect information to be returned
    params = {}
    params['teff'] = synth_model.teff()
    params['logg'] = synth_model.logg()
    params['MH'] = synth_model.MH()
    params['vmic'] = synth_model.vmic()
    params['vmac'] = synth_model.vmac()
    params['vsini'] = synth_model.vsini()
    params['limb_darkening_coeff'] = synth_model.limb_darkening_coeff()
    params['R'] = synth_model.R()
    params['continuum_correction'] = synth_model.continuum_correction()

    errors = {}
    errors['teff'] = synth_model.eteff()
    errors['logg'] = synth_model.elogg()
    errors['MH'] = synth_model.eMH()
    errors['vmic'] = synth_model.evmic()
    errors['vmac'] = synth_model.evmac()
    errors['vsini'] = synth_model.evsini()
    errors['limb_darkening_coeff'] = synth_model.elimb_darkening_coeff()
    errors['R'] = synth_model.eR()
    errors['continuum_correction'] = synth_model.econtinuum_correction()

    # Free abundances (original, transformed [X/H] [X/Fe] and errors)
    free_abundances = synth_model.transformed_free_abundances()

    status = {}
    status['days'] = synth_model.calculation_time.day-1
    status['hours'] = synth_model.calculation_time.hour
    status['minutes'] = synth_model.calculation_time.minute
    status['seconds'] = synth_model.calculation_time.second
    status['dof'] = synth_model.m.dof
    status['error'] = synth_model.m.errmsg
    status['chisq(tanh)'] = synth_model.chisq
    status['rchisq(tanh)'] = synth_model.reduced_chisq
    status['chisq'] = synth_model.basic_chisq
    status['rchisq'] = synth_model.reduced_basic_chisq
    status['niter'] = synth_model.m.niter
    status['nsynthesis'] = synth_model.m.nfev
    status['status'] = synth_model.m.status

    synth_spectrum = create_spectrum_structure(waveobs, synth_model.last_final_fluxes)
    if synth_model.last_continuum_correction is not None:
        # Uncorrect synthetic spectrum to recover the real theoretical spectrum
        synth_spectrum['flux'] *= synth_model.last_continuum_correction
        synth_spectrum['err'] *= synth_model.last_continuum_correction
        # Correct the observed spectrum
        spectrum['flux'] *= synth_model.last_continuum_correction
        spectrum['err'] *= synth_model.last_continuum_correction

    return spectrum, synth_spectrum, params, errors, free_abundances, status, stats_linemasks



class EquivalentWidthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, abundances, teff=5000, logg=3.0, MH=0.0, vmic=2.0):
        self.elements = {}
        #self.elements["1"] = "H"
        #self.elements["2"] = "He"
        self.elements["3"] = "Li"
        self.elements["4"] = "Be"
        self.elements["5"] = "B"
        self.elements["6"] = "C"
        self.elements["7"] = "N"
        self.elements["8"] = "O"
        self.elements["9"] = "F"
        self.elements["10"] = "Ne"
        self.elements["11"] = "Na"
        self.elements["12"] = "Mg"
        self.elements["13"] = "Al"
        self.elements["14"] = "Si"
        self.elements["15"] = "P"
        self.elements["16"] = "S"
        self.elements["17"] = "Cl"
        self.elements["18"] = "Ar"
        self.elements["19"] = "K"
        self.elements["20"] = "Ca"
        self.elements["21"] = "Sc"
        self.elements["22"] = "Ti"
        self.elements["23"] = "V"
        self.elements["24"] = "Cr"
        self.elements["25"] = "Mn"
        self.elements["26"] = "Fe"
        self.elements["27"] = "Co"
        self.elements["28"] = "Ni"
        self.elements["29"] = "Cu"
        self.elements["30"] = "Zn"
        self.elements["31"] = "Ga"
        self.elements["32"] = "Ge"
        self.elements["33"] = "As"
        self.elements["34"] = "Se"
        self.elements["35"] = "Br"
        self.elements["36"] = "Kr"
        self.elements["37"] = "Rb"
        self.elements["38"] = "Sr"
        self.elements["39"] = "Y"
        self.elements["40"] = "Zr"
        self.elements["41"] = "Nb"
        self.elements["42"] = "Mo"
        self.elements["43"] = "Tc"
        self.elements["44"] = "Ru"
        self.elements["45"] = "Rh"
        self.elements["46"] = "Pd"
        self.elements["47"] = "Ag"
        self.elements["48"] = "Cd"
        self.elements["49"] = "In"
        self.elements["50"] = "Sn"
        self.elements["51"] = "Sb"
        self.elements["52"] = "Te"
        self.elements["53"] = "I"
        self.elements["54"] = "Xe"
        self.elements["55"] = "Cs"
        self.elements["56"] = "Ba"
        self.elements["57"] = "La"
        self.elements["58"] = "Ce"
        self.elements["59"] = "Pr"
        self.elements["60"] = "Nd"
        self.elements["61"] = "Pm"
        self.elements["62"] = "Sm"
        self.elements["63"] = "Eu"
        self.elements["64"] = "Gd"
        self.elements["65"] = "Tb"
        self.elements["66"] = "Dy"
        self.elements["67"] = "Ho"
        self.elements["68"] = "Er"
        self.elements["69"] = "Tm"
        self.elements["70"] = "Yb"
        self.elements["71"] = "Lu"
        self.elements["72"] = "Hf"
        self.elements["73"] = "Ta"
        self.elements["74"] = "W"
        self.elements["75"] = "Re"
        self.elements["76"] = "Os"
        self.elements["77"] = "Ir"
        self.elements["78"] = "Pt"
        self.elements["79"] = "Au"
        self.elements["80"] = "Hg"
        self.elements["81"] = "Tl"
        self.elements["82"] = "Pb"
        self.elements["83"] = "Bi"
        self.elements["84"] = "Po"
        self.elements["85"] = "At"
        self.elements["86"] = "Rn"
        self.elements["87"] = "Fr"
        self.elements["88"] = "Ra"
        self.elements["89"] = "Ac"
        self.elements["90"] = "Th"
        self.elements["91"] = "Pa"
        self.elements["92"] = "U"
        self.elements["101"] = "Md"
        self.elements["106"] = "Sg"
        self.elements["107"] = "Bh"
        self.elements["108"] = "Hs"
        self.elements["112"] = "Cn"
        self.elements["113"] = "Uut"
        self.elements["114"] = "Uuq"

        self.modeled_layers_pack = modeled_layers_pack
        self.linelist = linelist
        self.abundances = abundances
        self.lines_for_teff = None
        self.lines_for_vmic = None
        #
        self.calculation_time = 0
        self.cache = {}
        p = [teff, logg, vmic]
        self._MH = MH
        self._eMH = 0.0

        self.min_MH = np.min(modeled_layers_pack[5])
        self.max_MH = np.max(modeled_layers_pack[5])
        #
        super(EquivalentWidthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p is not None:
            # Update internal structure for fitting:
            for i in xrange(len(p)):
                self._parinfo[i]['value'] = p[i]

        key = "%.2f %.2f %.2f %.2f " % (self.teff(), self.logg(), self.MH(), self.vmic())
        if self.cache.has_key(key):
            print "Cache:", key
            self.last_final_values = self.cache[key]
        else:
            print "Generating:", key
            # Optimization to avoid too small changes in parameters or repetition
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH())
            spec_abund, normal_abund, x_over_h, x_over_fe = determine_abundances(atmosphere_layers, \
                    self.teff(), self.logg(), self.MH(), self.linemasks, self.abundances, microturbulence_vel = self.vmic(), verbose=0)

            # First iteration
            if self.lines_for_teff is None or self.lines_for_vmic is None:
                self.select_good_lines(x_over_h, strict_teff=True, strict_vmic=True)

            values_to_evaluate = []
            fitted_lines_params = []
            selected_x_over_h = []
            for i, (lines_for_teff, lines_for_vmic) in enumerate(zip(self.lines_for_teff, self.lines_for_vmic)):
                ## Temperature
                # y = mx + c
                x = self.linemasks['lower state (eV)'][lines_for_teff]
                y = x_over_h[lines_for_teff]
                A = np.vstack([x, np.ones(len(x))]).T
                m1, c1 = np.linalg.lstsq(A, y)[0]
                #import matplotlib.pyplot as plt
                #plt.scatter(x, y)
                #plt.plot(x, m1*x + c1)
                #plt.show()

                ## Vmic
                # y = mx + c
                x = np.log10(self.linemasks['ew'][lines_for_vmic]/self.linemasks['wave_peak'][lines_for_vmic])
                y = x_over_h[lines_for_vmic]
                A = np.vstack([x, np.ones(len(x))]).T
                m2, c2 = np.linalg.lstsq(A, y)[0]
                #import matplotlib.pyplot as plt
                #plt.scatter(x, y)
                #plt.plot(x, m1*x + c1)
                #plt.show()


                ## Gravity
                abundance_diff = np.median(x_over_h[lines_for_teff]) - np.median(x_over_h[lines_for_vmic])
                #abundance_diff2 = self.MH() - np.median(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
                abundance_diff2 = self.MH() - np.median(x_over_h[self.lines_for_teff[0]]) # Always [0] where Fe 1 is

                print " # Element:                   ", self.teff_elements[i], self.vmic_elements[i]
                print "   Teff/Vmic slopes:            %.2f %.2f" % (m1, m2)
                print "   Abundances diff:             %.2f" % abundance_diff
                print "   Abundances diff with model:  %.2f" % abundance_diff2
                print "   Abundances stdev:            %.2f %.2f" % (np.std(x_over_h[lines_for_teff]), np.std(x_over_h[lines_for_vmic]))
                print "   Abundances median:           %.2f %.2f" % (np.median(x_over_h[lines_for_teff]), np.median(x_over_h[lines_for_vmic]))

                # Rounded to 3 and 2 decimals (using string convertion works better than np.round)
                values_to_evaluate.append(float("%.2f" % m1))
                values_to_evaluate.append(float("%.2f" % m2))
                values_to_evaluate.append(float("%.2f" % abundance_diff))
                #values_to_evaluate.append(float("%.2f" % abundance_diff2))
                #abundances_to_evaluate = np.arange(len(self.linemasks))
                ##abundances_to_evaluate[:] = 10.
                #abundances_to_evaluate[:] = 0.
                #abundances_to_evaluate[lines_for_teff] = x_over_h[lines_for_teff] - np.median(x_over_h[lines_for_teff])
                #abundances_to_evaluate[lines_for_vmic] = x_over_h[lines_for_vmic] - np.median(x_over_h[lines_for_vmic])
                # TODO: It will not work if there are more than 1 element (i.e. Fe and Ti)
                abundances_to_evaluate = x_over_h - np.median(x_over_h[np.logical_or(lines_for_teff, lines_for_vmic)])
                values_to_evaluate = np.hstack((values_to_evaluate,  abundances_to_evaluate)).tolist()
                #values_to_evaluate = abundances_to_evaluate

                residuals = np.asarray(values_to_evaluate) - self.y
                #print "Eval:", np.sum(np.tanh(self.weights*residuals)**2)
                print " - Chisq:                       %.4f" % np.sum((self.weights*residuals)**2)
                fitted_lines_params.append(m1)
                fitted_lines_params.append(c1)
                fitted_lines_params.append(m2)
                fitted_lines_params.append(c2)
                selected_x_over_h.append(lines_for_teff.copy())
                selected_x_over_h.append(lines_for_vmic.copy())
            self.last_final_values = (np.asarray(values_to_evaluate), x_over_h, selected_x_over_h, fitted_lines_params)
            self.cache[key] = self.last_final_values


        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        return values_to_evaluate.copy()

    # Default procedure to be called every iteration.  It simply prints
    # the parameter values.
    import scipy
    blas_enorm32, = scipy.lib.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float32))
    blas_enorm64, = scipy.lib.blas.get_blas_funcs(['nrm2'],np.array([0],dtype=np.float64))
    def defiter(self, fcn, x, iter, fnorm=None, functkw=None,
                       quiet=0, iterstop=None, parinfo=None,
                       format=None, pformat='%.10g', dof=1):

        if quiet:
            return
        if fnorm is None:
            [status, fvec] = fcn(x, fjac=None, **functkw)
            # If the returned fvec has more than four bits I assume that we have
            # double precision
            # It is important that the machar is determined by the precision of
            # the returned value, not by the precision of the input array
            if np.array([fvec]).dtype.itemsize>4:
                self.blas_enorm = mpfit.blas_enorm64
            else:
                self.blas_enorm = mpfit.blas_enorm32
            fnorm = self.enorm(fvec)**2

        # Determine which parameters to print
        nprint = len(x)
        print "*Iter ", ('%6i' % iter),"   CHI-SQUARE = ",('%.10g' % fnorm)," DOF = ", ('%i' % dof)
        for i in range(nprint):
            if (parinfo is not None) and (parinfo[i].has_key('parname')):
                p = '   ' + parinfo[i]['parname'] + ' = '
            else:
                p = '   P' + str(i) + ' = '
            if (parinfo is not None) and (parinfo[i].has_key('mpprint')):
                iprint = parinfo[i]['mpprint']
            else:
                iprint = 1
            if iprint:
                print p + (pformat % x[i]) + '  '
        ##### Lines
        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        self.select_good_lines(x_over_h, strict_teff=True, strict_vmic=True) # Modifies self.lines_for_teff and self.lines_for_vmic

        ##### Metallicity
        #self._MH = np.median(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        self._MH = np.median(x_over_h[self.lines_for_teff[0]]) # Only from Fe 1
        self._MH = np.min((self._MH, self.max_MH))
        self._MH = np.max((self._MH, self.min_MH))
        self._eMH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])

        return 0


    def select_good_lines(self, x_over_h, strict_teff=True, strict_vmic=False):
        """
            Modifies self.lines_for_teff and self.lines_for_vmic
        """
        # Out of range
        bad = np.logical_or(x_over_h > 1.0, x_over_h < -5)
        #### Line selection
        # Select elements for determining the Teff (traditionally Fe 1)
        self.lines_for_teff = []
        for element in self.teff_elements:
            lines_for_teff = self.linemasks['element'] == element

            if strict_teff and len(np.where(~bad & lines_for_teff)[0]) > 1:
                # Outliers
                x = self.linemasks['lower state (eV)'][~bad & lines_for_teff]
                y = x_over_h[~bad & lines_for_teff]
                A = np.vstack([x, np.ones(len(x))]).T
                m0, c0 = np.linalg.lstsq(A, y)[0]
                #lower_limit = m0*self.linemasks['lower state (eV)'] + c0 - 3*np.std(y)
                #upper_limit = m0*self.linemasks['lower state (eV)'] + c0 + 3*np.std(y)
                interq = np.percentile(y, 99) - np.percentile(y, 1)
                lower_limit = m0*self.linemasks['lower state (eV)'] + c0 - interq/2.
                upper_limit = m0*self.linemasks['lower state (eV)'] + c0 + interq/2.
                reject_lines_for_teff = np.logical_or(x_over_h > upper_limit, x_over_h < lower_limit)
                reject_lines_for_teff = np.logical_or(reject_lines_for_teff, bad)
                #import matplotlib.pyplot as plt
                #plt.scatter(x, y)
                #plt.plot(x, m0*x + c0)
                #plt.plot(x, m0*x + c0 + 3*np.std(y))
                #plt.plot(x, m0*x + c0 - 3*np.std(y))
                #plt.show()

                # Discard bad lines and outliers
                clean_lines_for_teff = np.logical_and(~reject_lines_for_teff, lines_for_teff)
                # Ensure that there are at least some lines
                if len(np.where(clean_lines_for_teff)[0]) <= 1:
                    # Discard only bad lines
                    clean_lines_for_teff = lines_for_teff[~bad]
                    if len(np.where(clean_lines_for_teff)[0]) <= 1:
                        clean_lines_for_teff = lines_for_teff
                if len(np.where(clean_lines_for_teff)[0]) <= 1:
                    raise Exception("Not enought lines for Teff (%i lines)" % len(np.where(clean_lines_for_teff)[0]))
                else:
                    self.lines_for_teff.append(clean_lines_for_teff)
            else:
                ##### ACCEPT all
                if len(np.where(~bad & lines_for_teff)[0]) > 0:
                    self.lines_for_teff.append(np.logical_and(lines_for_teff, np.logical_not(bad)))
                else:
                    self.lines_for_teff.append(lines_for_teff)
            print " > Selected", element, "lines for teff:", len(np.where(self.lines_for_teff[-1])[0]), "of", len(np.where(lines_for_teff)[0])


        # Select elements for determining the Vmic (traditionally Fe 2)
        self.lines_for_vmic = []
        for element in self.vmic_elements:
            lines_for_vmic = self.linemasks['element'] == element


            if strict_vmic and len(np.where(~bad & lines_for_vmic)[0]) > 1:
                x = np.log10(self.linemasks['ew'][~bad & lines_for_vmic]/self.linemasks['wave_peak'][~bad & lines_for_vmic])
                y = x_over_h[~bad & lines_for_vmic]
                A = np.vstack([x, np.ones(len(x))]).T
                m0, c0 = np.linalg.lstsq(A, y)[0]
                #lower_limit = m0*np.log10(self.linemasks['ew']/self.linemasks['wave_peak']) + c0 - 3*np.std(y)
                #upper_limit = m0*np.log10(self.linemasks['ew']/self.linemasks['wave_peak']) + c0 + 3*np.std(y)
                interq = np.percentile(y, 99) - np.percentile(y, 1)
                lower_limit = m0*np.log10(self.linemasks['ew']/self.linemasks['wave_peak']) + c0 - interq/2.
                upper_limit = m0*np.log10(self.linemasks['ew']/self.linemasks['wave_peak']) + c0 + interq/2.
                reject_lines_for_vmic = np.logical_or(x_over_h > upper_limit, x_over_h < lower_limit)
                reject_lines_for_vmic = np.logical_or(reject_lines_for_vmic, bad)
                #import matplotlib.pyplot as plt
                #plt.scatter(x, y)
                #plt.plot(x, m0*x + c0)
                #plt.plot(x, m0*x + c0 + 3*np.std(y))
                #plt.plot(x, m0*x + c0 - 3*np.std(y))
                #plt.show()

                #discard = np.logical_or(bad, reject1)
                #discard = np.logical_or(discard, reject2)

                # Discard bad lines and outliers
                clean_lines_for_vmic = np.logical_and(~reject_lines_for_vmic, lines_for_vmic)
                if len(np.where(clean_lines_for_vmic)[0]) <= 1:
                    clean_lines_for_vmic = lines_for_vmic[~bad]
                    # Discard only bad lines
                    if len(np.where(clean_lines_for_vmic)[0]) <= 1:
                        clean_lines_for_vmic = lines_for_vmic
                if len(np.where(clean_lines_for_vmic)[0]) <= 1:
                    raise Exception("Not enought lines for Vmic calculations (%i lines)" % len(np.where(clean_lines_for_vmic)[0]))
                else:
                    self.lines_for_vmic.append(clean_lines_for_vmic)
            else:
                ##### ACCEPT all
                if len(np.where(~bad & lines_for_vmic)[0]) > 0:
                    self.lines_for_vmic.append(np.logical_and(lines_for_vmic, np.logical_not(bad)))
                else:
                    self.lines_for_vmic.append(lines_for_vmic)
            print " > Selected", element, "lines for vmic:", len(np.where(self.lines_for_vmic[-1])[0]), "of", len(np.where(lines_for_vmic)[0])



    def fitData(self, linemasks, teff_elements=["Fe 1"], vmic_elements=["Fe 2"], parinfo=None, max_iterations=20, quiet=True):
        base = 3
        if len(parinfo) < base:
            raise Exception("Wrong number of parameters!")

        if len(teff_elements) != len(vmic_elements):
            raise Exception("Inconsistent number of teff/vmic elements!")

        if teff_elements[0] != "Fe 1" or vmic_elements[0] != "Fe 2":
            raise Exception("First element should be always Fe 1/2!")


        if sys.platform == "win32":
            # On Windows, the best timer is time.clock()
            default_timer = time.clock
        else:
            # On most other platforms the best timer is time.time()
            default_timer = time.time
        self.linemasks = linemasks
        self.teff_elements = teff_elements
        self.vmic_elements = vmic_elements
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 0.0   # Residuals are limited between -1.0 and 1.0 (np.tanh(residuals/1.0)) minimizing the influence of bad fluxes
                     # * Spectrum must be a normalized one
        #chisq_limit = 0.0002 # np.sum(np.asarray([0.00, 0.00, 0.01, 0.01])**2))
        chisq_limit = 0.0001 # np.sum(np.asarray([0.00, 0.00, 0.01])**2))
        _t0 = default_timer()

        #index = np.asarray([0, 1, 2])
        #index = np.arange(3*len(teff_elements)) # 3 values: zero slopes and zero difference between element1 and element2
        index = np.arange(3*len(teff_elements) + len(linemasks)) # 3 values: zero slopes and zero difference between element1 and element2
        #index = np.arange(len(linemasks))
        #index = np.arange(4*len(teff_elements)) # 3 values: zero slopes and zero difference between element1 and element2
        target_values = np.zeros(len(index))
        weights = np.ones(len(index))
        super(EquivalentWidthModel, self).fitData(index, target_values, weights=weights, parinfo=parinfo, chisq_limit=chisq_limit, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=max_iterations, quiet=quiet, iterfunct=self.defiter)

        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        residuals = values_to_evaluate - target_values
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
        # Chisq using tanh
        self.chisq = np.sum(np.tanh(weights * residuals)**2)
        self.reduced_chisq = self.chisq / self.m.dof
        # Chisq without using tanh for minimizing outliers
        self.basic_chisq = np.sum((weights * residuals)**2)
        self.reduced_basic_chisq = self.basic_chisq / self.m.dof

        self.cache = {}

        _t1 = default_timer()
        sec = timedelta(seconds=int(_t1 - _t0))
        self.calculation_time = datetime(1,1,1) + sec

    def teff(self): return self._parinfo[0]['value']
    def logg(self): return self._parinfo[1]['value']
    def vmic(self): return self._parinfo[2]['value']

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def evmic(self): return self.m.perror[2]

    def MH(self): return self._parinfo[3]['value']
    def eMH(self): return self.m.perror[3]

    #def MH(self): return self._MH
    #def eMH(self): return self._eMH

    def print_solution(self):
        # Calculate MH
        values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = self.last_final_values
        #MH = np.median(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        #MH = np.std(x_over_h[np.logical_or(self.lines_for_teff[0], self.lines_for_vmic[0])])
        MH = np.median(x_over_h[self.lines_for_teff[0]]) # Only from Fe 1
        eMH = np.std(x_over_h[self.lines_for_teff[0]]) # Only from Fe 1

        header = "%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","vmic")
        #solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), self.MH(), self.vmic())
        #errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), self.eMH(), self.evmic())
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.teff(), self.logg(), MH, self.vmic())
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f" % (self.eteff(), self.elogg(), eMH, self.evmic())

        print "           ", header
        print "Solution:  ", solution
        print "Errors:    ", errors
        print ""

        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","chisq(tanh)","rchisq(tanh)","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.4f\t%8.4f\t%8.4f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.chisq, self.reduced_chisq, self.basic_chisq, self.reduced_basic_chisq, self.rms)
        print ""
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status


def modelize_spectrum_from_EW(linemasks, modeled_layers_pack, linelist, abundances, initial_teff, initial_logg, initial_MH, initial_vmic, teff_elements=["Fe 1"], vmic_elements=["Fe 2"], max_iterations=20):
    """
    """
    teff_range = modeled_layers_pack[3]
    logg_range = modeled_layers_pack[4]
    MH_range = modeled_layers_pack[5]

    parinfo = __create_EW_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, teff_range, logg_range, MH_range)

    EW_model = EquivalentWidthModel(modeled_layers_pack, linelist, abundances, MH=initial_MH)

    lfilter = None
    for element in np.hstack((teff_elements, vmic_elements)):
        if lfilter is None:
            lfilter = linemasks['element'] == element
        else:
            lfilter = np.logical_or(lfilter, linemasks['element'] == element)
    linemasks = linemasks[lfilter]
    EW_model.fitData(linemasks, teff_elements=teff_elements, vmic_elements=vmic_elements, parinfo=parinfo, max_iterations=max_iterations, quiet=False)
    print "\n"
    EW_model.print_solution()

    # Collect information to be returned
    params = {}
    params['teff'] = EW_model.teff()
    params['logg'] = EW_model.logg()
    params['MH'] = EW_model.MH()
    params['vmic'] = EW_model.vmic()

    errors = {}
    errors['teff'] = EW_model.eteff()
    errors['logg'] = EW_model.elogg()
    errors['MH'] = EW_model.eMH()
    errors['vmic'] = EW_model.evmic()

    status = {}
    status['days'] = EW_model.calculation_time.day-1
    status['hours'] = EW_model.calculation_time.hour
    status['minutes'] = EW_model.calculation_time.minute
    status['seconds'] = EW_model.calculation_time.second
    status['dof'] = EW_model.m.dof
    status['error'] = EW_model.m.errmsg
    status['chisq(tanh)'] = EW_model.chisq
    status['rchisq(tanh)'] = EW_model.reduced_chisq
    status['chisq'] = EW_model.basic_chisq
    status['rchisq'] = EW_model.reduced_basic_chisq
    status['niter'] = EW_model.m.niter
    status['nsynthesis'] = EW_model.m.nfev
    status['status'] = EW_model.m.status

    values_to_evaluate, x_over_h, selected_x_over_h, fitted_lines_params = EW_model.last_final_values

    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params

