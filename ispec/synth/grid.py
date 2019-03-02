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
import time
import numpy as np
from scipy import spatial
import pandas as pd
import multiprocessing
from multiprocessing import Pool
from lockfile import FileLock, LockTimeout, AlreadyLocked
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, Column
import cPickle as pickle
import logging

from ispec.atmospheres import interpolate_atmosphere_layers, valid_atmosphere_target, _interpolate
from ispec.common import mkdir_p, estimate_vmic, estimate_vmac
from ispec.spectrum import create_spectrum_structure, resample_spectrum, create_wavelength_filter, read_spectrum, write_spectrum
from effects import apply_post_fundamental_effects
from ispec.modeling.common import Constants
import ispec.synth.common

def generate_fundamental_spectrum(grid, waveobs, teff, logg, MH, alpha, microturbulence_vel, regions=None):
    """
    Generates an interpolated spectrum from a grid. vmic is always fixed and depends
    on the grid existing points.

    No macroturbulence, rotation (vsini), limb darkening coefficient or resolution is considered
    in this process. That's why it is named as "fundamental" spectrum.
    """
    return generate_spectrum(grid, waveobs, teff, logg, MH, alpha, microturbulence_vel, macroturbulence=0.0, vsini=0.0, limb_darkening_coeff=0.00, R=0, regions=regions)

def _add_target_if_possible(free_parameters, target_point, name, value):
    if name in free_parameters:
        target_point.append(value)
    return target_point

def generate_spectrum(grid, waveobs, teff, logg, MH, alpha, microturbulence_vel, macroturbulence=0.0, vsini=0.0, limb_darkening_coeff=0.00, R=0, regions=None):
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulation, kdtree, ranges, base_dirname = grid

    if regions is None:
        global_wave_base = np.min(waveobs)
        global_wave_top = np.max(waveobs)
        regions = np.recarray((1,),  dtype=[('wave_base', float), ('wave_top', float)])
        regions['wave_base'][0] = global_wave_base
        regions['wave_top'][0] = global_wave_top
    # Only read the part of the spectrum that needs to be interpolated
    def read_point_value(f, regions=None):
        return read_spectrum(f, apply_filters=False, sort=False, regions=regions)
    custom_read_point_value = lambda f: read_point_value(f, regions=regions)

    target_point = []
    target_point = _add_target_if_possible(free_parameters, target_point, 'teff', teff)
    target_point = _add_target_if_possible(free_parameters, target_point, 'logg', logg)
    target_point = _add_target_if_possible(free_parameters, target_point, 'MH', MH)
    target_point = _add_target_if_possible(free_parameters, target_point, 'alpha', alpha)
    target_point = _add_target_if_possible(free_parameters, target_point, 'vmic', microturbulence_vel)
    interpolated = _interpolate(delaunay_triangulation, kdtree, existing_points, filenames, custom_read_point_value, value_fields, target_point)
    interpolated_spectrum = create_spectrum_structure(interpolated['waveobs'], interpolated['flux'])

    # Make sure we return the number of expected fluxes
    if not np.array_equal(interpolated_spectrum['waveobs'], waveobs):
        interpolated_spectrum = resample_spectrum(interpolated_spectrum, waveobs, method="linear", zero_edges=True)

    segments = None
    vrad = (0,)
    interpolated_spectrum['flux'] = apply_post_fundamental_effects(interpolated_spectrum['waveobs'], interpolated_spectrum['flux'], segments, \
                    macroturbulence=macroturbulence, vsini=vsini, \
                    limb_darkening_coeff=limb_darkening_coeff, R=R, vrad=vrad)
    return interpolated_spectrum['flux']


def __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, resolution, pickled_modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, code="spectrum", use_molecules=False, tmp_dir=None, locked=False):
    multiprocessing.current_process().daemon=False

    import dill # To allow pickle of lambda functions (e.g., one element in modeled_layers_pack)
    import pickle
    modeled_layers_pack = pickle.loads(pickled_modeled_layers_pack)

    if valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
        if not locked:
            lock = FileLock(filename_out+".lock")
            try:
                lock.acquire(timeout=-1)    # Don't wait
            except (LockTimeout, AlreadyLocked) as e:
                # Some other process is computing this spectrum, do not continue
                print "Skipping", teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, "already locked"
                return None

        try:
            print "[started]", teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, resolution
            # Prepare atmosphere model
            atmosphere_layers = interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)
            fixed_abundances=None
            # Synthesis
            synth_spectrum = create_spectrum_structure(wavelengths)
            synth_spectrum['flux'] = ispec.synth.common.generate_spectrum(synth_spectrum['waveobs'], \
                    atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
                    fixed_abundances, microturbulence_vel = vmic, \
                    macroturbulence=vmac, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
                    R=resolution, regions=segments, verbose=0, \
                    code=code, use_molecules=use_molecules, tmp_dir=tmp_dir)
            # FITS
            write_spectrum(synth_spectrum, filename_out)
            print "[finished]", teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, resolution
        finally:
            if not locked: # Not locked in this function
                lock.release()
    else:
        raise Exception("Not valid: %i %.2f %.2f" % (teff, logg, MH))



def precompute_synthetic_grid(output_dirname, ranges, wavelengths, to_resolution, modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, segments=None, number_of_processes=1, code="spectrum", use_molecules=False, steps=False, tmp_dir=None):
    """
    Pre-compute a synthetic grid with some reference ranges (Teff, log(g) and
    MH combinations) and all the steps that iSpec will perform in the
    astrophysical parameter determination process.

    All the non-convolved spectra will be saved in a subdir and a complete
    grid file with the reference points already convolved will be saved in a
    FITS file for fast comparison.

    The output directory can be used by the routines 'model_spectrum' and
    'estimate_initial_ap'.
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'synthe', 'sme']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    reference_list_filename = output_dirname + "/parameters.tsv"
    if to_resolution is not None:
        reference_grid_filename = output_dirname + "/convolved_grid_%i.fits.gz" % to_resolution
    fits_dir = os.path.join(output_dirname, "grid/")
    mkdir_p(fits_dir)
    if steps:
        steps_fits_dir = os.path.join(output_dirname, "steps/")
        mkdir_p(steps_fits_dir)

    import dill # To allow pickle of lambda functions (e.g., one element in modeled_layers_pack)
    import pickle
    pickled_modeled_layers_pack = pickle.dumps(modeled_layers_pack)

    # For code != "grid", ranges are always in position 7 (for grid it would be in position 8)
    valid_ranges = modeled_layers_pack[7]
    teff_range = valid_ranges['teff']
    logg_range = valid_ranges['logg']
    MH_range = valid_ranges['MH']
    alpha_range = valid_ranges.get('alpha', (-1.5, 1.5)) # Default (0.,) if 'alpha' is not a free parameter for atmosphere interpolation
    vmic_range = valid_ranges.get('vmic', (0.0, 50.)) # Default (0.,) if 'vmic' is not a free parameter for atmosphere interpolation

    # Parallelization pool
    if number_of_processes == 1:
        pool = None
    else:
        pool = Pool(number_of_processes)

    # Create grid binary file
    elapsed = 0 # seconds

    num_ref_spec = len(ranges)
    num_spec = num_ref_spec * 9 # Reference + 8 variations in Teff, logg, MH, alpha, vmic, vmac, vsini, limb darkening coeff

    i = 0
    for teff, logg, MH, alpha, vmic in ranges:
        if vmic is None:
            vmic = estimate_vmic(teff, logg, MH)
        vmac = 0.0 # This can be modified after synthesis if needed
        vsini = 0.0 # This can be modified after synthesis if needed
        limb_darkening_coeff = 0.00 # This can be modified after synthesis if needed
        resolution = 0 # This can be modified after synthesis if needed
        is_step = False
        if not valid_atmosphere_target(modeled_layers_pack, {'teff': teff, 'logg': logg, 'MH': MH, 'alpha': alpha}):
            raise Exception("Target parameters out of the valid ranges: teff={} logg={} MH={} alpha={}".format(teff, logg, MH, alpha))
        points = [
                    (teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step),
                ]
        if steps:
            is_step = True
            new_teff = teff+Constants.SYNTH_STEP_TEFF if teff+Constants.SYNTH_STEP_TEFF <= teff_range[-1] else teff-Constants.SYNTH_STEP_TEFF
            new_logg = logg+Constants.SYNTH_STEP_LOGG if logg+Constants.SYNTH_STEP_LOGG <= logg_range[-1] else logg-Constants.SYNTH_STEP_LOGG
            new_MH = MH+Constants.SYNTH_STEP_MH if MH+Constants.SYNTH_STEP_MH <= MH_range[-1] else MH-Constants.SYNTH_STEP_MH
            new_alpha = alpha+Constants.SYNTH_STEP_ALPHA if alpha+Constants.SYNTH_STEP_ALPHA <= alpha_range[-1] else alpha-Constants.SYNTH_STEP_ALPHA
            new_vmic = vmic+Constants.SYNTH_STEP_VMIC if vmic+Constants.SYNTH_STEP_VMIC <= vmic_range[-1] else vmic-Constants.SYNTH_STEP_VMIC
            # For each reference point, calculate also the variations that iSpec will perform in the first iteration
            points += [ # Final unconvolved spectra where vmic/vmac are free and do not follow vmic/vmac empirical relations
                        (new_teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, new_logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, logg, new_MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, logg, MH, new_alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, logg, MH, alpha, new_vmic, vmac, vsini, limb_darkening_coeff, is_step),
                    ]
            points += [
                        # Final unconvolved spectra where vmic is not free and does follow vmic empirical relations
                        (new_teff, logg, MH, alpha, estimate_vmic(new_teff, logg, MH), vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, new_logg, MH, alpha, estimate_vmic(teff, new_logg, MH), vmac, vsini, limb_darkening_coeff, is_step),
                        (teff, logg, new_MH, alpha, estimate_vmic(teff, logg, new_MH), vmac, vsini, limb_darkening_coeff, is_step),
                    ]

        for j, (teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, is_step) in enumerate(points):
            if is_step:
                filename_out = steps_fits_dir + "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}".format(int(teff), logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff) + ".fits.gz"
            else:
                filename_out = fits_dir + "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}".format(int(teff), logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff) + ".fits.gz"

            if os.path.exists(filename_out):
                print "Skipping", teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, "already computed"
                continue


            if pool is None:
                if sys.platform == "win32":
                    # On Windows, the best timer is time.clock()
                    default_timer = time.clock
                else:
                    # On most other platforms the best timer is time.time()
                    default_timer = time.time

                lock = FileLock(filename_out+".lock")
                try:
                    lock.acquire(timeout=-1)    # Don't wait
                except (LockTimeout, AlreadyLocked) as e:
                    # Some other process is computing this spectrum, do not continue
                    print "Skipping", teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, "already locked"
                    continue

                try:
                    tcheck = default_timer()
                    # Validate parameters
                    __generate_synthetic_fits(filename_out, wavelengths, segments, teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, resolution, pickled_modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, code=code, use_molecules=use_molecules, tmp_dir=tmp_dir, locked=True)
                    elapsed = default_timer() - tcheck

                    print "-----------------------------------------------------"
                    print "Remaining time:"
                    print "\t", (num_spec-i)*elapsed, "seconds"
                    print "\t", (num_spec-i)*(elapsed/60), "minutes"
                    print "\t", (num_spec-i)*(elapsed/(60*60)), "hours"
                    print "\t", (num_spec-i)*(elapsed/(60*60*24)), "days"
                    print "-----------------------------------------------------"
                finally:
                    lock.release()

            else:
                pool.apply_async(__generate_synthetic_fits, [filename_out, wavelengths, segments, teff, logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff, resolution, pickled_modeled_layers_pack, atomic_linelist, isotopes, solar_abundances], kwds={'code': code, 'use_molecules': use_molecules, 'tmp_dir':tmp_dir, 'locked':False})
            i += 1

    if pool is not None:
        pool.close()
        pool.join()



    # Create parameters.tsv
    reference_list = Table()
    if len(np.unique(ranges[['logg', 'MH', 'alpha', 'vmic']])) == len(ranges):
        reference_list.add_column(Column(name='fixed_teff', dtype=int))
    else:
        reference_list.add_column(Column(name='teff', dtype=int))
    if len(np.unique(ranges[['teff', 'MH', 'alpha', 'vmic']])) == len(ranges):
        reference_list.add_column(Column(name='fixed_logg', dtype=float))
    else:
        reference_list.add_column(Column(name='logg', dtype=float))
    if len(np.unique(ranges[['teff', 'logg', 'alpha', 'vmic']])) == len(ranges):
        reference_list.add_column(Column(name='fixed_MH', dtype=float))
    else:
        reference_list.add_column(Column(name='MH', dtype=float))
    if len(np.unique(ranges[['teff', 'logg', 'MH', 'vmic']])) == len(ranges):
        reference_list.add_column(Column(name='fixed_alpha', dtype=float))
    else:
        reference_list.add_column(Column(name='alpha', dtype=float))
    if len(np.unique(ranges[['teff', 'logg', 'MH', 'alpha']])) == len(ranges):
        reference_list.add_column(Column(name='fixed_vmic', dtype=float))
    else:
        reference_list.add_column(Column(name='vmic', dtype=float))
    reference_list.add_column(Column(name='filename', dtype='|S100'))
    for teff, logg, MH, alpha, vmic in ranges:
        # Only use the first spectra generated for each combination
        zero_vmac = 0.0
        zero_vsini = 0.0
        zero_limb_darkening_coeff = 0.00
        reference_filename_out = "./grid/{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}".format(int(teff), logg, MH, alpha, vmic, zero_vmac, zero_vsini, zero_limb_darkening_coeff) + ".fits.gz"
        reference_list.add_row((int(teff), logg, MH, alpha, vmic, reference_filename_out))

    if not os.path.exists(reference_list_filename):
        lock = FileLock(reference_list_filename+".lock")
        try:
            lock.acquire(timeout=-1)    # Don't wait
        except (LockTimeout, AlreadyLocked) as e:
            # Some other process is writing this file, do not continue
            print "Skipping", reference_list_filename, "already locked"
        else:
            try:
                ascii.write(reference_list, reference_list_filename, delimiter='\t', overwrite=True)
                print "Written", reference_list_filename
            finally:
                lock.release()

    if to_resolution is not None:
        if not os.path.exists(reference_grid_filename):
            lock = FileLock(reference_grid_filename+".lock")
            try:
                lock.acquire(timeout=-1)    # Don't wait
            except (LockTimeout, AlreadyLocked) as e:
                # Some other process is computing this spectrum, do not continue
                print "Skipping", reference_grid_filename, "already locked"
            else:
                try:
                    reference_grid = None
                    complete_reference_list = Table()
                    complete_reference_list.add_column(Column(name='teff', dtype=int))
                    complete_reference_list.add_column(Column(name='logg', dtype=float))
                    complete_reference_list.add_column(Column(name='MH', dtype=float))
                    complete_reference_list.add_column(Column(name='alpha', dtype=float))
                    complete_reference_list.add_column(Column(name='vmic', dtype=float))
                    complete_reference_list.add_column(Column(name='vmac', dtype=float))
                    complete_reference_list.add_column(Column(name='vsini', dtype=float))
                    complete_reference_list.add_column(Column(name='limb_darkening_coeff', dtype=float))
                    for teff, logg, MH, alpha, vmic in ranges:
                        # Only use the first spectra generated for each combination
                        zero_vmac = 0.0
                        zero_vsini = 0.0
                        zero_limb_darkening_coeff = 0.00
                        vmac = estimate_vmac(teff, logg, MH)
                        vsini = 1.6 # Sun
                        limb_darkening_coeff = 0.6
                        reference_filename_out = "{0}_{1:.2f}_{2:.2f}_{3:.2f}_{4:.2f}_{5:.2f}_{6:.2f}_{7:.2f}".format(int(teff), logg, MH, alpha, vmic, zero_vmac, zero_vsini, zero_limb_darkening_coeff) + ".fits.gz"
                        if not os.path.exists(fits_dir + reference_filename_out):
                            continue
                        complete_reference_list.add_row((int(teff), logg, MH, alpha, vmic, vmac, vsini, limb_darkening_coeff))


                        # Spectra in the grid is convolved to the specified resolution for fast comparison
                        print "Quick grid:", reference_filename_out
                        spectrum = read_spectrum(fits_dir + reference_filename_out)

                        segments = None
                        vrad = (0,)
                        spectrum['flux'] = apply_post_fundamental_effects(spectrum['waveobs'], spectrum['flux'], segments, \
                                    macroturbulence=vmac, vsini=vsini, \
                                    limb_darkening_coeff=limb_darkening_coeff, R=to_resolution, vrad=vrad)

                        if reference_grid is None:
                            reference_grid = spectrum['flux']
                        else:
                            reference_grid = np.vstack((reference_grid, spectrum['flux']))

                    if len(ranges) == len(complete_reference_list):
                        # Generate FITS file with grid for fast comparison
                        primary_hdu = fits.PrimaryHDU(reference_grid)
                        wavelengths_hdu = fits.ImageHDU(wavelengths, name="WAVELENGTHS")
                        params_bintable_hdu = fits.BinTableHDU(complete_reference_list.as_array(), name="PARAMS")
                        fits_format = fits.HDUList([primary_hdu, wavelengths_hdu, params_bintable_hdu])
                        fits_format.writeto(reference_grid_filename, overwrite=True)
                        print "Written", reference_grid_filename
                finally:
                    lock.release()

def estimate_initial_ap(spectrum, precomputed_dir, resolution, linemasks, default_teff = 5000., default_logg = 2.5, default_MH = 0.0, default_alpha = 0.0, default_vmic = 1.0, default_vmac = 0.0, default_vsini = 0.0, default_limb_darkening_coeff = 0.00):
    """
    Estimate the initial atmospheric parameters by using a pre-computed grid
    at a given resolution. The comparison will be based on the linemasks.
    """

    reference_grid_filename = precomputed_dir + "/convolved_grid_%i.fits.gz" % resolution
    estimation_found = False
    if not os.path.exists(reference_grid_filename):
        logging.warn("Pre-computed grid does not exists for R = %i" % resolution)
    else:
        try:
            grid = fits.open(reference_grid_filename)
            grid_waveobs = np.asarray(grid['WAVELENGTHS'].data, dtype=float)
            resampled_spectrum = resample_spectrum(spectrum, grid_waveobs, method="linear")

            fsegment = create_wavelength_filter(resampled_spectrum, regions=linemasks)
            fsegment = np.logical_and(fsegment, resampled_spectrum['flux'] > 0.0)
            isegment = np.where(fsegment == True)[0]

            # http://en.wikipedia.org/wiki/Goodness_of_fit#Example
            residuals = grid['PRIMARY'].data[:,isegment] - resampled_spectrum['flux'][isegment]
            chisq = np.sum((residuals)**2, axis=1)
            min_j = np.argmin(chisq)
            initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff = grid['PARAMS'].data[min_j]
            estimation_found = True

        except Exception, e:
            print "Initial parameters could not be estimated"
            print type(e), e.message
            pass
        finally:
            grid.close()

    if estimation_found:
        return initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff
    else:
        return default_teff, default_logg, default_MH, default_alpha, default_vmic, default_vmac, default_vsini, default_limb_darkening_coeff

def load_spectral_grid(input_path):
    """
    :param input_path:
        Name of the input file (i.e. models.dump) or directory (new interpolator)
    :type input_path: string

    :returns:
        List of modeled_layers, used_values_for_layers, proximity, teff_range, logg_range, MH_range and nlayers
    """
    if not os.path.isdir(input_path):
        raise Exception("Input path '{}' is not a directory".format(input_path))

    base_dirname = input_path
    atm_dirname = os.path.join(base_dirname, "grid")
    params_filename = os.path.join(base_dirname, "parameters.tsv")
    cache_filename = os.path.join(base_dirname, "cache.dump")

    if not os.path.exists(atm_dirname):
        raise Exception("Grid path '{}' does not exist".format(atm_dirname))
    if not os.path.exists(params_filename):
        raise Exception("Parameters file '{}' does not exist".format(params_filename))

    parameters = pd.read_csv(params_filename, sep="\t")
    # Divide grid in two using the temperature to reduce the number of data points
    # to be used by a single Delaunay Triangulation
    teff_limit = 10000
    teff_margin = teff_limit/2
    filter_cool_parameters = np.array(parameters['teff'] < teff_limit)
    filter_intermediate_parameters = np.array(np.logical_and(parameters['teff'] > teff_limit-teff_margin, parameters['teff'] < teff_limit+teff_margin))
    filter_hot_parameters = np.array(parameters['teff'] >= teff_limit)
    parameters_subsets = [filter_cool_parameters, filter_intermediate_parameters, filter_hot_parameters]
    #
    filenames = base_dirname + "/" + parameters['filename']
    filenames = np.asarray(filenames)
    ranges = {}
    free_parameters = parameters.columns.get_values().tolist()
    free_parameters = filter(lambda x: not x.startswith("fixed_") and x != "filename", free_parameters)
    for free_param in free_parameters:
        free_param_range = np.unique(parameters[free_param])
        ranges[free_param] = free_param_range
    if len(free_parameters) == 0:
        raise Exception("No free parameters in grid: '{}'".format(input_path))
    existing_points = parameters[free_parameters]
    existing_points = np.array(existing_points)

    # The delaunay triangulation and kdtree can be computationally expensive,
    # do it once and save a pickle dump
    # But first verify that the computed one (if exists) matches what is expected
    use_dump = False
    if os.path.exists(cache_filename):
        use_dump = True
        delaunay_triangulations, kdtree = pickle.load(open(cache_filename, 'rb'))
        for delaunay_triangulation, parameters_subset in zip(delaunay_triangulations['precomputed'], parameters_subsets):
            if delaunay_triangulation is None:
                continue
            if not (delaunay_triangulation.ndim == len(free_parameters) and delaunay_triangulation.npoints == len(existing_points[parameters_subset])):
                use_dump = False
                break

    if not os.path.exists(cache_filename) or not use_dump:
        delaunay_triangulations = {'subsets': parameters_subsets, 'precomputed': []}
        for i, parameters_subset in enumerate(parameters_subsets):
            logging.info("Pre-computing [{}/{}]...".format(i+1, len(parameters_subsets)))
            if len(existing_points[parameters_subset]) > 0:
                delaunay_triangulations['precomputed'].append(spatial.Delaunay(existing_points[parameters_subset]))
            else:
                delaunay_triangulations['precomputed'].append(None)
        kdtree = spatial.cKDTree(existing_points)
        pickle.dump((delaunay_triangulations, kdtree), open(cache_filename, 'wb'), protocol=2)

    # Functions will receive the parameters in the same order
    #read_point_value = lambda f, regions: read_spectrum(f, apply_filters=False, sort=False, regions=regions)
    def read_point_value(f, regions=None):
        return read_spectrum(f, apply_filters=False, sort=False, regions=regions)

    value_fields = ["waveobs", "flux"]
    return existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname

def valid_interpolated_spectrum_target(grid, target):
    """
        Returns False if model could not be interpolated
    """
    existing_points, free_parameters, filenames, read_point_value, value_fields, delaunay_triangulations, kdtree, ranges, base_dirname = grid
    target_point = []
    for param in free_parameters:
        target_point.append(target[param])
    target_point_cannot_be_interpolated = None
    for delaunay_triangulation in delaunay_triangulations['precomputed']:
        if delaunay_triangulation is None:
            continue
        simplex = delaunay_triangulation.find_simplex(target_point)
        if target_point_cannot_be_interpolated is None:
            target_point_cannot_be_interpolated = np.any(simplex == -1)
        else:
            target_point_cannot_be_interpolated = target_point_cannot_be_interpolated and np.any(simplex == -1)
    return not target_point_cannot_be_interpolated if target_point_cannot_be_interpolated is not None else False

