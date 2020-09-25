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
import os
cimport numpy as np


cdef extern from "stdio.h":
    ctypedef struct FILE

cdef extern from "spectrum/spectrum.h":
    struct memo:
        int lyman
        int balmer
        int paschen
        int brackett
        int pfund
        int humphreys
        int hprofl
        int helium
        int strong
        int interval
    ctypedef struct linedata:
        double wave
        double code
        int iso
        double atomass
        double abund
        double chi1
        double chi2
        double chi3
        double chi4
        double chi
        double Eu
        double El
        double gf
        double wavel
        double waveh
        float  xnum[100]
        float  a[100]
        float  dopp[100]
        float  capnu[100]
        float  dlg[100]
        char   T[5]
        double alp
        double sig
        double gammar
        double gammas
        double gammaw
        double fac
        int ai
        int flag

cdef extern from "synthesizer_func.h":
    ctypedef void (*progressfunc)(double num, void *user_data)
    int ew_and_depth(char *atmosphere_model_file, char *linelist_file, char *isotope_file, char *abundances_file, double microturbulence_vel, double start, double end, int verbose, int num_measures, double *output_wave, double *output_code, double *output_ew, double *output_depth, progressfunc user_func, void *user_data)
    int synthesize_spectrum(char *atmosphere_model_file, char *linelist_file, char *isotope_file, char *abundances_file, char *fixed_abundances_file, double microturbulence_vel, int verbose, int num_measures, double* waveobs, double* waveobs_mask, double *fluxes, progressfunc user_func, void *user_data)
    int macroturbulence_spectrum(double *waveobs, double *fluxes, int num_measures, double macroturbulence, int verbose, progressfunc user_func, void *user_data)
    int rotation_spectrum(double *waveobs, double *fluxes, int num_measures, double vsini, double limb_darkening_coeff, int verbose, progressfunc user_func, void *user_data)
    int resolution_spectrum(double *waveobs, double *fluxes, int num_measures, int R, int verbose, progressfunc user_func, void *user_data)
    int abundances_determination(char *atmosphere_model_file, char *linelist_file, int num_measures, char *abundances_file, double microturbulence_vel, int verbose, double* ignore, double* abundances, double *normal_abundances, double*relative_abundances, progressfunc user_func, void *user_data)

#### Callback
def dummy_func(double num):
    pass

cdef void callback(double num, void *f):
    (<object>f)(num)
##############

# waveobs in armstrong
# microtturbulence velocity in km/s
def spectrum(np.ndarray[np.double_t,ndim=1] waveobs, np.ndarray[np.double_t,ndim=1] waveobs_mask, char* atmosphere_model_file, char* linelist_file = "input/linelists/default.300_1100nm.lst", char *isotope_file = "input/abundances/isotope.iso", char* abundances_file = "input/abundances/default.stdatom.dat", char* fixed_abundances_file="none", double microturbulence_vel = 2.0, double macroturbulence = 3.0, double vsini = 2.0, double limb_darkening_coeff = 0.0, int R=500000, int nlayers = 56, int verbose = 0, update_progress_func=None):
    if not os.path.exists(atmosphere_model_file):
        raise Exception("Atmosphere model file '%s' does not exists!" % atmosphere_model_file)
    if not os.path.exists(linelist_file):
        raise Exception("Line list file '%s' does not exists!" % linelist_file)
    if not os.path.exists(abundances_file):
        raise Exception("Abundances file '%s' does not exists!" % abundances_file)
    global Ntau
    global flagr
    global flagc
    global flagk
    global flagg
    global flagmph
    global flagI
    global flagt
    global flagp
    global flagP
    global flagu
    global flagO
    global flagC
    global mghla
    global mghlb
    global mu
    global NI
    Ntau = nlayers  # 72 layers for castelli-kurucz atmosphere models, 56 for MARCS
    flagr = 0
    flagc = 0
    flagk = 0
    flagg = 0
    flagmgh = 0
    flagI = 1   # Isotopes (1: True, 0: False)
    flagt = 0
    flagp = 0
    flagP = 0
    flagu = 0
    flagO = 0
    flagC = 0
    mghla = 0
    mghlb = 0
    mu = 1.0
    NI = 0

    cdef int num_measures = len(waveobs)
    cdef np.ndarray[np.double_t,ndim=1] fluxes = np.zeros(num_measures, dtype=float)
    if num_measures <= 1:
        # We need at least 2 wavelengths, if not return an zeroed result
        return fluxes

    if update_progress_func==None:
        update_progress_func = dummy_func

    synthesize_spectrum(atmosphere_model_file, linelist_file, isotope_file, abundances_file,
            fixed_abundances_file,
            microturbulence_vel, verbose, num_measures, <double*> waveobs.data,
            <double*> waveobs_mask.data,
            <double*> fluxes.data, callback, <void*>update_progress_func)

    fluxes = apply_post_fundamental_effects(waveobs, fluxes, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, R, verbose, update_progress_func)
    return fluxes



def apply_post_fundamental_effects(np.ndarray[np.double_t,ndim=1] waveobs, np.ndarray[np.double_t,ndim=1] fluxes, double microturbulence_vel = 2.0, double macroturbulence = 3.0, double vsini = 2.0, double limb_darkening_coeff = 0.0, int R=500000, int verbose = 0, update_progress_func=None):

    cdef int num_measures = len(waveobs)
    if num_measures <= 1:
        # We need at least 2 wavelengths, if not return the same fluxes
        return fluxes

    if update_progress_func==None:
        update_progress_func = dummy_func

    if macroturbulence > 0:
        macroturbulence_spectrum(<double*> waveobs.data, <double*> fluxes.data,
            num_measures, macroturbulence, verbose, callback, <void*>update_progress_func)
    if vsini > 0 or limb_darkening_coeff > 0:
        rotation_spectrum(<double*> waveobs.data, <double*> fluxes.data, num_measures,
            vsini, limb_darkening_coeff, verbose, callback, <void*>update_progress_func)
    if R > 0:
        resolution_spectrum(<double*> waveobs.data, <double*> fluxes.data,
            num_measures, R, verbose, callback, <void*>update_progress_func)
    return fluxes



# microtturbulence velocity in km/s
def abundances(char* atmosphere_model_file, char* linelist_file, int num_measures, np.ndarray[np.double_t,ndim=1] ignore, char* abundances_file, double microturbulence_vel = 2.0, int nlayers=56, int verbose = 0, update_progress_func=None):
    if not os.path.exists(atmosphere_model_file):
        raise Exception("Atmosphere model file '%s' does not exists!" % atmosphere_model_file)
    if not os.path.exists(linelist_file):
        raise Exception("Line list file '%s' does not exists!" % linelist_file)
    if not os.path.exists(abundances_file):
        raise Exception("Abundances file '%s' does not exists!" % abundances_file)
    global Ntau
    global flagr
    global flagc
    global flagk
    global flagg
    global flagmph
    global flagI
    global flagt
    global flagp
    global flagP
    global flagu
    global flagO
    global flagC
    global mghla
    global mghlb
    global mu
    global NI
    global flagCNO
    Ntau = nlayers  # 72 layers for castelli-kurucz atmosphere models, 56 for MARCS
    flagr = 0
    flagc = 0
    flagk = 0
    flagg = 0
    flagmgh = 0
    flagI = 1   # Isotopes (1: True, 0: False)
    flagt = 0
    flagp = 0
    flagP = 0
    flagu = 0
    flagO = 0
    flagC = 0
    mghla = 0
    mghlb = 0
    mu = 1.0
    NI = 0
    flagCNO = 0

    cdef np.ndarray[np.double_t,ndim=1] abundances = np.zeros(num_measures, dtype=float)
    cdef np.ndarray[np.double_t,ndim=1] normal_abundances = np.zeros(num_measures, dtype=float)
    cdef np.ndarray[np.double_t,ndim=1] relative_abundances = np.zeros(num_measures, dtype=float)
    if num_measures <= 0:
        # We need at least 2 wavelengths, if not return an zeroed result
        return abundances, normal_abundances, relative_abundances

    if update_progress_func==None:
        update_progress_func = dummy_func

    abundances_determination(atmosphere_model_file, linelist_file, num_measures, abundances_file,
            microturbulence_vel, verbose,
            <double*> ignore.data,
            <double*> abundances.data,
            <double*> normal_abundances.data, <double*> relative_abundances.data,
            callback, <void*>update_progress_func)

    return abundances, normal_abundances, relative_abundances


def calculate_ew_and_depth(char* atmosphere_model_file, char* linelist_file, char *isotope_file, char* abundances_file, int num_lines, double microturbulence_vel = 2.0, int nlayers = 56, double start=3000, double end=11000, int verbose = 0, update_progress_func=None):
    if not os.path.exists(atmosphere_model_file):
        raise Exception("Atmosphere model file '%s' does not exists!" % atmosphere_model_file)
    if not os.path.exists(linelist_file):
        raise Exception("Line list file '%s' does not exists!" % linelist_file)
    if not os.path.exists(abundances_file):
        raise Exception("Abundances file '%s' does not exists!" % abundances_file)
    global Ntau
    global flagr
    global flagc
    global flagk
    global flagg
    global flagmph
    global flagI
    global flagt
    global flagp
    global flagP
    global flagu
    global flagO
    global flagC
    global mghla
    global mghlb
    global mu
    global NI
    Ntau = nlayers  # 72 layers for castelli-kurucz atmosphere models, 56 for MARCS
    flagr = 0
    flagc = 0
    flagk = 0
    flagg = 0
    flagmgh = 0
    flagI = 1   # Isotopes (1: True, 0: False)
    flagt = 0
    flagp = 0
    flagP = 0
    flagu = 0
    flagO = 0
    flagC = 0
    mghla = 0
    mghlb = 0
    mu = 1.0
    NI = 0

    cdef np.ndarray[np.double_t,ndim=1] output_wave = np.zeros(num_lines, dtype=float)
    cdef np.ndarray[np.double_t,ndim=1] output_code = np.zeros(num_lines, dtype=float)
    cdef np.ndarray[np.double_t,ndim=1] output_ew = np.zeros(num_lines, dtype=float)
    cdef np.ndarray[np.double_t,ndim=1] output_depth= np.zeros(num_lines, dtype=float)
    if num_lines <= 0:
        # We need at least 1 wavelengths, if not return an zeroed result
        return output_wave, output_ew, output_depth

    if update_progress_func==None:
        update_progress_func = dummy_func

    ew_and_depth(atmosphere_model_file, linelist_file, isotope_file, abundances_file, \
            microturbulence_vel, start, end, verbose, num_lines, \
            <double*> output_wave.data, <double*> output_code.data, <double*> output_ew.data, <double*> output_depth.data, \
           callback, <void*>update_progress_func)

    return output_wave, output_code, output_ew, output_depth

