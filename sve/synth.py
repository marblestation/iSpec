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
import os
import sys
import time
from datetime import datetime, timedelta
import numpy as np
from mpfitmodels import *
from abundances import *
from atmospheres import *
from spectrum import *
from multiprocessing import Process
from multiprocessing import Queue
from Queue import Empty

import log
import logging

def read_SPECTRUM_linelist(linelist_filename):
    """
    Load a SPECTRUM linelist for spectral synthesis
    """
    linelist = np.array([tuple(line.rstrip('\r\n').split()) for line in open(linelist_filename,)], \
                            dtype=[('wave (A)', '<f8'), ('species', '|S10'), ('lower state (cm^-1)', int), \
                            ('upper state (cm^-1)', int), ('log(gf)', '<f8'), ('fudge factor', '<f8'), \
                            ('transition type', '|S10'), ('rad', '<f8'),  ('stark', '<f8'), ('waals', '<f8'), \
                            ('note', '|S100')])
    return linelist

def write_SPECTRUM_linelist(linelist, linelist_filename=None):
    """
    Saves a SPECTRUM linelist for spectral synthesis.
    If filename is not specified, a temporary file is created and the name is returned.
    """
    if linelist_filename != None:
        out = open(linelist_filename, "w")
    else:
        # Temporary file
        out = tempfile.NamedTemporaryFile(delete=False)
    out.write("\n".join(["  ".join(map(str, (line['wave (A)'], line['species'], line['lower state (cm^-1)'], line['upper state (cm^-1)'], line['log(gf)'], line['fudge factor'], line['transition type'], line['rad'], line['stark'], line['waals'], line['note']))) for line in linelist]))
    out.close()
    return out.name


def generate_fundamental_spectrum(waveobs, waveobs_mask, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel, verbose=0, update_progress_func=None, timeout=600, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs only
    if waveobs_mask contains the value 1.0 at the same position.

    No macroturbulence, rotation (vsini), limb darkening coefficient or resolution is considered
    in this process. That's why it is named as "fundamental" spectrum.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be empty.
    """
    return generate_spectrum(waveobs, waveobs_mask, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence = 0.0, vsini = 0.0, limb_darkening_coeff = 0.0, R=0, verbose=verbose, update_progress_func=update_progress_func, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file)


def generate_spectrum(waveobs, waveobs_mask, atmosphere_layers, teff, logg, MH, linelist, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, update_progress_func=None, timeout=600, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None):
    """
    Generates a synthetic spectrum for the wavelength specified in waveobs only
    if waveobs_mask contains the value 1.0 at the same position.

    The atmosphere model, linelist, abundances and fixed abundances can be specified
    as numpy recarray tables and they will be saved to disk to be used by SPECTRUM. In case the
    user already has the information saved onto the disk, the filenames can be
    specified to reduce input/output time (and the numpy recarray tables will be ignored)

    Fixed abundances can be empty.
    """

    # OPTIMIZATION: Use already saved files to reduce input/output time
    remove_tmp_atm_file = False
    remove_tmp_abund_file = False
    remove_tmp_fixed_abund_file = False
    remove_tmp_linelist_file = False
    if atmosphere_layers_file == None:
        atmosphere_layers_file = write_atmosphere(atmosphere_layers, teff, logg, MH)
        remove_tmp_atm_file = True
    if abundances_file == None:
        abundances_file = write_SPECTRUM_abundances(abundances)
        remove_tmp_abund_file = True
    if fixed_abundances_file == None:
        fixed_abundances_file = write_SPECTRUM_fixed_abundances(fixed_abundances)
        remove_tmp_fixed_abund_file = True
    if linelist_file == None:
        linelist_file = write_SPECTRUM_linelist(linelist)
        remove_tmp_linelist_file = True
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

    p = Process(target=__generate_spectrum, args=(result_queue, waveobs, waveobs_mask, atmosphere_layers_file, linelist_file, abundances_file, fixed_abundances_file), kwargs={'microturbulence_vel': microturbulence_vel, 'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'nlayers': nlayers, 'verbose': verbose, 'update_progress_func':update_progress_func})
    p.start()
    fluxes = np.zeros(len(waveobs))
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            fluxes = result_queue.get(timeout=1)
        except Empty:
            # No results, continue waiting
            num_seconds += 1
        else:
            # Results received!
            break
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

def __generate_spectrum(result_queue, waveobs, waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, nlayers=56, verbose=0, update_progress_func=None):
    """
    Generate synthetic spectrum and apply macroturbulence, rotation (visini), limb darkening coeff and resolution except
    if all those parameters are set to zero, in that case the fundamental synthetic spectrum is returned.
    """
    import synthesizer
    fluxes = synthesizer.spectrum(waveobs*10., waveobs_mask, atmosphere_model_file, linelist_file, abundances_file, fixed_abundances_file, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, nlayers, verbose, update_progress_func)
    result_queue.put(fluxes)



def apply_post_fundamental_effects(waveobs, fluxes, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, update_progress_func=None, timeout=600):
    """
    Apply macroturbulence, rotation (vsini), limb darkening coefficient and/or resolution
    """

    # Generate spectrum should be run in a separate process in order
    # to force the reload of the "synthesizer" module which
    # contains C code with static variables in functions that should
    # be reinitialized to work properly
    # * The best solution would be to improve the C code but since it is too complex
    #   this hack has been implemented
    result_queue = Queue()

    # TODO: Allow communications between process in order to update the GUI progress bar
    update_progress_func = None

    p = Process(target=__apply_post_fundamental_effects, args=(result_queue, waveobs, fluxes), kwargs={'macroturbulence': macroturbulence, 'vsini': vsini, 'limb_darkening_coeff': limb_darkening_coeff, 'R': R, 'verbose': verbose, 'update_progress_func':update_progress_func})
    p.start()
    num_seconds = 0
    # Constantly check that the process has not died without returning any result and blocking the queue call
    while p.is_alive() and num_seconds < timeout:
        try:
            fluxes = result_queue.get(timeout=1)
        except Empty:
            # No results, continue waiting
            num_seconds += 1
        else:
            # Results received!
            break
    if num_seconds >= timeout:
        logging.error("A timeout has occurred in the application of post fundamental effects.")
        p.terminate()
    elif np.all(fluxes == 0):
        logging.error("The application of post fundamental effects has failed.")
        p.terminate()
    else:
        p.join()

    return fluxes


def __apply_post_fundamental_effects(result_queue, waveobs, fluxes, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.0, R=500000, verbose=0, update_progress_func=None):
    """
    Apply macroturbulence, rotation (visini), limb darkening coeff and resolution to already generated fundamental synthetic spectrum.
    """
    import synthesizer
    fluxes = synthesizer.apply_post_fundamental_effects(waveobs*10., fluxes, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, verbose, update_progress_func)
    result_queue.put(fluxes)


class SynthModel(MPFitModel):
    """
    Match synthetic spectrum to observed spectrum
    * Requires the synthetic spectrum generation functionality on
    """
    def __init__(self, modeled_layers_pack, linelist, abundances, fixed_abundances, teff=5000, logg=3.0, MH=0.0, vmic=2.0, vmac=0.0, vsini=2.0, limb_darkening_coeff=0.0, R=0):
        self.modeled_layers_pack = modeled_layers_pack
        self.linelist = linelist
        self.abundances = abundances
        self.fixed_abundances = fixed_abundances
        #
        self.calculation_time = 0
        self.waveobs = None
        self.waveobs_mask = None
        self.cache = {}
        p = [teff, logg, MH, vmic, vmac, vsini, limb_darkening_coeff, R]
        #
        self.abundances_file = None
        self.linelist_file = None
        super(SynthModel, self).__init__(p)

    def _model_function(self, x, p=None):
        # The model function with parameters p required by mpfit library
        if p != None:
            # Update internal structure for fitting:
            self._parinfo[0]['value'] = p[0]
            self._parinfo[1]['value'] = p[1]
            self._parinfo[2]['value'] = p[2]
            self._parinfo[3]['value'] = p[3]
            self._parinfo[4]['value'] = p[4]
            self._parinfo[5]['value'] = p[5]
            self._parinfo[6]['value'] = p[6]
            self._parinfo[7]['value'] = p[7]

        complete_key = "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %i" % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()))
        key = "%.2f %.2f %.2f %.2f" % (self.teff(), self.logg(), self.MH(), self.vmic())
        if self.cache.has_key(key):
            print "Cache:", complete_key
            self.last_fluxes = self.cache[key].copy()
        else:
            print "Generating:", complete_key
            # No fixed abundances
            fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])
            # Atmosphere
            atmosphere_layers = interpolate_atmosphere_layers(self.modeled_layers_pack, self.teff(), self.logg(), self.MH())
            # Fundamental synthetic fluxes
            self.last_fluxes = generate_fundamental_spectrum(self.waveobs, self.waveobs_mask, atmosphere_layers, self.teff(), self.logg(), self.MH(), self.linelist, self.abundances, self.fixed_abundances, microturbulence_vel=self.vmic(), abundances_file=self.abundances_file, linelist_file=self.linelist_file, verbose=0)
            # Optimization to avoid too small changes in parameters or repetition
            self.cache[key] = self.last_fluxes.copy()

        self.last_final_fluxes = apply_post_fundamental_effects(self.waveobs, self.last_fluxes, macroturbulence=self.vmac(), vsini=self.vsini(), limb_darkening_coeff=self.limb_darkening_coeff(), R=0, verbose=0)
        if self.R()>0:
            synth_spectrum = np.recarray((len(self.waveobs), ), dtype=[('waveobs', float),('flux', float),('err', float)])
            synth_spectrum['waveobs'] = self.waveobs
            synth_spectrum['flux'] = self.last_final_fluxes
            synth_spectrum['err'] = 0.0
            convolved_synth_spectrum = convolve_spectrum(synth_spectrum, self.R())
            self.last_final_fluxes = convolved_synth_spectrum['flux']
        return self.last_final_fluxes[self.comparing_mask]

    def fitData(self, waveobs, waveobs_mask, comparing_mask, fluxes, weights=None, parinfo=None, quiet=True):
        if len(self._parinfo) != 8:
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
        if weights == None:
            weights = np.ones(len(waveobs))
        ftol = 1.e-4 # Terminate when the improvement in chisq between iterations is ftol > -(new_chisq/chisq)**2 +1
        xtol = 1.e-4
        gtol = 1.e-4
        damp = 1.0   # Residuals are limited between -1.0 and 1.0 (np.tanh(residuals/1.0)) minimizing the influence of bad fluxes
                     # * Spectrum must be a normalized one
        maxiter = 20 # Maximum number of iterations
        _t0 = default_timer()

        # Write abundances and linelist to avoid writing the same info in each iteration
        self.abundances_file = write_SPECTRUM_abundances(self.abundances)
        self.linelist_file = write_SPECTRUM_linelist(self.linelist)

        super(SynthModel, self).fitData(waveobs[self.comparing_mask], fluxes[self.comparing_mask], weights=weights[self.comparing_mask], parinfo=parinfo, ftol=ftol, xtol=xtol, gtol=gtol, damp=damp, maxiter=maxiter, quiet=quiet)

        residuals = self.last_final_fluxes[self.comparing_mask] - fluxes[self.comparing_mask]
        self.rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
        # Chisq without using tanh
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

    def eteff(self): return self.m.perror[0]
    def elogg(self): return self.m.perror[1]
    def eMH(self): return self.m.perror[2]
    def evmic(self): return self.m.perror[3]
    def evmac(self): return self.m.perror[4]
    def evsini(self): return self.m.perror[5]
    def elimb_darkening_coeff(self): return self.m.perror[6]
    def eR(self): return self.m.perror[7]

    def print_solution(self):
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("teff","logg","MH","vmic","vmac","vsini","limb","R")
        solution = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.teff(), self.logg(), self.MH(), self.vmic(), self.vmac(), self.vsini(), self.limb_darkening_coeff(), int(self.R()))
        errors = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8i" % (self.eteff(), self.elogg(), self.eMH(), self.evmic(), self.evmac(), self.evsini(), self.elimb_darkening_coeff(), int(self.eR()))
        print "         ", header
        print "Solution:", solution
        print "Errors:  ", errors
        print "Calculation time:\t%d:%d:%d:%d" % (self.calculation_time.day-1, self.calculation_time.hour, self.calculation_time.minute, self.calculation_time.second)
        header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("DOF","niter","nsynthesis","chisq(tanh)","rchisq(tanh)","chisq","rchisq","rms")
        stats = "%8i\t%8i\t%8i\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (self.m.dof, self.m.niter, self.m.nfev, self.chisq, self.reduced_chisq, self.basic_chisq, self.reduced_basic_chisq, self.rms)
        print "         ", header
        print "Stats:   ", stats
        print "Return code:", self.m.status

def __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, teff_range, logg_range, MH_range):
    """
    Creates the structure needed for the mpfitmodel
    """
    free_params = [param.lower() for param in free_params]
    parinfo = [{'value':0., 'fixed':False, 'limited':[False, False], 'limits':[0., 0.], 'step':0} for i in np.arange(8)]
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
    parinfo[1]['step'] = 0.50 # For auto-derivatives
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
    parinfo[5]['step'] = 0.5 # For auto-derivatives
    parinfo[5]['limited'] = [True, True]
    parinfo[5]['limits'] = [0.0, 50.0]
    #
    parinfo[6]['parname'] = "limb_darkening_coef"
    parinfo[6]['value'] = initial_limb_darkening_coeff
    parinfo[6]['fixed'] = not parinfo[6]['parname'].lower() in free_params
    parinfo[6]['step'] = 0.10 # For auto-derivatives
    parinfo[6]['limited'] = [True, True]
    parinfo[6]['limits'] = [0.0, 1.0]
    #
    parinfo[7]['parname'] = "R"
    parinfo[7]['value'] = initial_R
    parinfo[7]['fixed'] = not parinfo[7]['parname'].lower() in free_params
    parinfo[7]['step'] = 100.0 # For auto-derivatives
    parinfo[7]['limited'] = [True, True]
    parinfo[7]['limits'] = [0.0, 300000.0]

    return parinfo


def __filter_linelist(linelist, segments):
    # Build wavelength points from regions
    lfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if lfilter == None:
            lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
        else:
            lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.))

    if lfilter != None:
        return linelist[lfilter]
    else:
        return linelist

def __create_waveobs_mask(waveobs, segments):
    # Build wavelength points from regions
    wfilter = None
    for region in segments:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if wfilter == None:
            wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
        else:
            wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
    waveobs_mask = np.zeros(len(waveobs))
    waveobs_mask[wfilter] = 1.0 # Compute fluxes only for selected segments

    return waveobs_mask

def __create_comparing_mask(waveobs, linemasks):
    # Build wavelength points from regions
    wfilter = None
    for region in linemasks:
        wave_base = region['wave_base']
        wave_top = region['wave_top']

        if wfilter == None:
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

def modelize_spectrum(spectrum, modeled_layers_pack, linelist, abundances, fixed_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, segments=None, linemasks=None):

    waveobs = spectrum['waveobs']
    flux = spectrum['flux']
    #if np.all(np.logical_and(spectrum['err'] > 0.0, spectrum['err'] <= 1.0)):
        #weights = 1./spectrum['err']
    #else:
        #weights = np.ones(len(waveobs))
    weights = np.ones(len(waveobs))
    #weights = 1. / (flux * 0.10)
    #snr = 100
    #weights = 1. / (flux / snr)
    #weights = 1. / spectrum['err']

    if segments == None:
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

    if linemasks == None:
        comparing_mask = np.ones(len(waveobs)) # Compare all fluxes
    else:
        comparing_mask = __create_comparing_mask(waveobs, linemasks)

    teff_range = modeled_layers_pack[3]
    logg_range = modeled_layers_pack[4]
    MH_range = modeled_layers_pack[5]

    parinfo = __create_param_structure(initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, free_params, teff_range, logg_range, MH_range)

    synth_model = SynthModel(modeled_layers_pack, linelist, abundances, fixed_abundances)
    synth_model.fitData(waveobs, waveobs_mask, comparing_mask, flux, weights=weights, parinfo=parinfo, quiet=False)
    stats_linemasks = __get_stats_per_linemask(waveobs, flux, synth_model.last_final_fluxes, weights, free_params, linemasks, verbose=True)
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

    errors = {}
    errors['teff'] = synth_model.eteff()
    errors['logg'] = synth_model.elogg()
    errors['MH'] = synth_model.eMH()
    errors['vmic'] = synth_model.evmic()
    errors['vmac'] = synth_model.evmac()
    errors['vsini'] = synth_model.evsini()
    errors['limb_darkening_coeff'] = synth_model.elimb_darkening_coeff()
    errors['R'] = synth_model.eR()

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

    synth_spectrum = np.recarray((len(waveobs), ), dtype=[('waveobs', float),('flux', float),('err', float)])
    synth_spectrum['waveobs'] = waveobs
    synth_spectrum['flux'] = synth_model.last_final_fluxes
    synth_spectrum['err'] = 0.0

    return synth_spectrum, params, errors, status, stats_linemasks

