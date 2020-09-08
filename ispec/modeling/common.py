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
import numpy as np



class Constants(object):
    ###################################
    # CONSTANTS
    ###################################
    SYNTH_STEP_TEFF = 100.
    SYNTH_STEP_LOGG = 0.1
    SYNTH_STEP_MH = 0.05
    SYNTH_STEP_ALPHA = 0.05
    SYNTH_STEP_VMIC = 0.5
    SYNTH_STEP_VMAC = 2.0
    SYNTH_STEP_VSINI = 2.0
    SYNTH_STEP_LIMB_DARKENING_COEFF = 0.05
    SYNTH_STEP_R = 100
    SYNTH_STEP_VRAD = 2
    SYNTH_STEP_ABUNDANCES = 0.05
    SYNTH_STEP_LOGGF = 0.01
    ###################################
    EW_STEP_TEFF = 500.
    EW_STEP_LOGG = 0.5
    EW_STEP_MH = 0.05
    EW_STEP_ALPHA = 0.05
    EW_STEP_VMIC = 0.5

def _filter_linemasks_not_in_segments(linemasks, segments):
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

def _create_comparing_mask(waveobs, linemasks, segments):
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

def _get_stats_per_linemask(waveobs, fluxes, synthetic_fluxes, weights, free_params, linemasks, verbose=False):

    results = np.recarray((len(linemasks), ), dtype=[('wave_peak', float),('wave_base', float),('wave_top', float),('chisq', float),('rchisq', float),('wchisq', float),('rwchisq', float),('rms', float)])
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

            # Unweighted
            chisq = np.sum((residuals)**2)
            reduced_chisq = chisq / dof

            # Weighted
            wchisq = np.sum((weights[wfilter] * residuals)**2)
            reduced_wchisq = wchisq / dof
        else:
            rms = -9999
            chisq = -9999
            reduced_chisq = -9999
            wchisq = -9999
            reduced_wchisq = -9999

        results['rms'][i] = rms
        # Unweighted
        results['chisq'][i] = chisq
        results['rchisq'][i] = reduced_chisq
        # Weighted
        results['wchisq'][i] = wchisq
        results['rwchisq'][i] = reduced_wchisq
        if verbose:
            header = "%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s" % ("wave_peak","wave_base","wave_top","wchisq","rwchisq","chisq","rchisq","rms")
            stats = "%8.2f\t%8.2f\t%8.2f\t%8.2f\t%8.4f\t%8.2f\t%8.4f\t%8.4f" % (wave_peak, wave_base, wave_top, wchisq, reduced_wchisq, chisq, reduced_chisq, rms)
            if i == 0:
                print("         ", header)
            print("Line     ", stats)
        i += 1

    return results

