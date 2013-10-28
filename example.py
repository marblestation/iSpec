#!/usr/bin/env python
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
import logging

################################################################################
#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/sblancoc/shared/iSpec/'):
    # avakas
    ispec_dir = '/home/sblancoc/shared/iSpec/'
elif os.path.exists('/home/blanco/shared/iSpec/'):
    # vanoise
    ispec_dir = '/home/blanco/shared/iSpec/'
else:
    ispec_dir = '/home/marble/shared/iSpec/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


#--- Change LOG level ----------------------------------------------------------
#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))
################################################################################


def read_write_spectrum():
    #--- Reading spectra -----------------------------------------------------------
    logging.info("Reading spectra")
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    ispec.write_spectrum(sun_spectrum, "example_sun.s")
    return sun_spectrum


def convert_air_to_vacuum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Converting wavelengths from air to vacuum and viceversa -------------------
    sun_spectrum_vacuum = ispec.air_to_vacuum(sun_spectrum)
    sun_spectrum_air = ispec.vacuum_to_air(sun_spectrum_vacuum)
    return sun_spectrum_vacuum, sun_spectrum_air

def plot():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Plotting (requires graphical interface) -----------------------------------
    logging.info("Plotting...")
    ispec.plot_spectra([sun_spectrum, mu_cas_a_spectrum])
    ispec.show_histogram(sun_spectrum['flux'])

def cut_spectrum_from_range():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Cut -----------------------------------------------------------------------
    logging.info("Cutting...")

    # - Keep points between two given wavelengths
    wfilter = ispec.create_wavelength_filter(sun_spectrum, wave_base=480.0, wave_top=670.0)
    cutted_sun_spectrum = sun_spectrum[wfilter]
    return cutted_sun_spectrum

def cut_spectrum_from_segments():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Cut -----------------------------------------------------------------------
    logging.info("Cutting...")
    # Keep only points inside a list of segments
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    wfilter = ispec.create_wavelength_filter(sun_spectrum, regions=segments)
    cutted_sun_spectrum = sun_spectrum[wfilter]
    return cutted_sun_spectrum


def determine_radial_velocity_with_mask():
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Radial Velocity determination with linelist mask --------------------------
    logging.info("Radial velocity determination with linelist mask...")
    # - Read atomic data
    mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Arcturus.372_926nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Sun.372_926nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.G2.375_679nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K0.378_679nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/Synthetic.Sun.300_1100nm.txt"
    #mask_file = ispec_dir + "input/linelists/CCF/VALD.Sun.300_1100nm.txt"
    ccf_mask = ispec.read_linelist_mask(mask_file)

    xcoord, fluxes, errors = ispec.build_velocity_profile(mu_cas_a_spectrum, \
                                                linelist=ccf_mask, lower_velocity_limit=-200.0, \
                                                upper_velocity_limit=200.0, velocity_step=1.0)

    models = ispec.modelize_velocity_profile(xcoord, fluxes, errors)
    best = ispec.select_good_velocity_profile_models(models, xcoord, fluxes)
    models = models[best]

    # Number of models represent the number of components
    components = len(models)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    return rv, components


def determine_radial_velocity_with_template():
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Radial Velocity determination with template -------------------------------
    logging.info("Radial velocity determination with template...")
    # - Read synthetic template
    template = ispec.read_spectrum(ispec_dir + \
            "/input/spectra/synthetic/Synth_ATLAS9.APOGEE_VALD_5777.0_4.44_0.0_1.0.txt.gz")
    # - Read observed template
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    xcoord, fluxes, errors = ispec.build_velocity_profile(mu_cas_a_spectrum, \
                                                template=template, lower_velocity_limit=-200.0, \
                                                upper_velocity_limit=200.0, velocity_step=1.0)

    models = ispec.modelize_velocity_profile(xcoord, fluxes, errors)
    best = ispec.select_good_velocity_profile_models(models, xcoord, fluxes)
    models = models[best]

    # Number of models represent the number of components
    components = len(models)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    return rv, components

def correct_radial_velocity():
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Radial Velocity correction ------------------------------------------------
    logging.info("Radial velocity correction...")
    rv = -96.40 # km/s
    mu_cas_a_spectrum = ispec.correct_velocity(mu_cas_a_spectrum, rv)


def determine_tellurics_shift():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Telluric velocity shift determination from spectrum --------------------------
    logging.info("Telluric velocity shift determination...")
    # - Telluric
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"
    linelist_telluric = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    xcoord, fluxes, errors = ispec.build_velocity_profile(sun_spectrum, \
                                            linelist=linelist_telluric, lower_velocity_limit=-100.0, \
                                            upper_velocity_limit=100.0, velocity_step=0.5)

    models = ispec.modelize_velocity_profile(xcoord, fluxes, errors, only_one_peak=True)
    bv = np.round(models[0].mu(), 2) # km/s
    return bv


def degrade_resolution():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Resolution degradation ----------------------------------------------------
    logging.info("Resolution degradation...")
    from_resolution = 80000
    to_resolution = 40000
    convolved_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, to_resolution, \
                                                    from_resolution=from_resolution)
    return convolved_sun_spectrum


def resample_spectrum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Resampling  --------------------------------------------------------------
    logging.info("Resampling...")
    wavelengths = np.arange(480.0, 670.0, 0.001)
    resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths)
    return resampled_sun_spectrum

def coadd_spectra():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Resampling and combining --------------------------------------------------
    logging.info("Resampling and comibining...")
    wavelengths = np.arange(480.0, 670.0, 0.001)
    resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths)
    resampled_mu_cas_a_spectrum = ispec.resample_spectrum(mu_cas_a_spectrum, wavelengths)
    # Coadd previously resampled spectra
    coadded_spectrum = ispec.create_spectrum_structure(resampled_sun_spectrum['waveobs'])
    coadded_spectrum['flux'] = resampled_sun_spectrum['flux'] + resampled_mu_cas_a_spectrum['flux']
    coadded_spectrum['err'] = np.sqrt(np.power(resampled_sun_spectrum['err'],2) + \
                                    np.power(resampled_mu_cas_a_spectrum['err'],2))
    return coadded_spectrum



def normalize_spectrum_using_continuum_regions():
    """
    Consider only continuum regions for the fit, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    continuum_regions = ispec.read_continuum_regions(ispec_dir + "/input/regions/fe_lines_continuum.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                            continuum_regions=continuum_regions, nknots=nknots, degree=degree, \
                            median_wave_range=median_wave_range, \
                            max_wave_range=max_wave_range, \
                            model=model, order=order, \
                            automatic_strong_line_detection=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_spectrum_in_segments():
    """
    Fit continuum in each segment independently, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = 1
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                            independent_regions=segments, nknots=nknots, degree=degree,\
                            median_wave_range=median_wave_range, \
                            max_wave_range=max_wave_range, \
                            model=model, order=order, \
                            automatic_strong_line_detection=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy2():
    """
    Use the whole spectrum, strategy 'max+median'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first MAXIMUM values and secondly medians in order to find the continuum
    order='max+median'
    median_wave_range=3.0
    max_wave_range=0.5

    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy1():
    """
    Use the whole spectrum, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy1_ignoring_strong_lines():
    """
    Use the whole spectrum but ignoring some strong lines, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")

    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def filter_cosmic_rays():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Filtering cosmic rays -----------------------------------------------------
    # Spectrum should be already normalized
    cosmics = ispec.create_filter_cosmic_rays(sun_spectrum, sun_continuum_model, \
                                            resampling_wave_step=0.001, window_size=15, \
                                            variation_limit=0.01)
    clean_sun_spectrum = sun_spectrum[~cosmics]
    return clean_sun_spectrum


def find_continuum_regions():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Find continuum regions ----------------------------------------------------
    logging.info("Finding continuum regions...")
    resolution = 80000
    sigma = 0.001
    max_continuum_diff = 1.0
    fixed_wave_step = 0.05
    sun_continuum_regions = ispec.find_continuum(sun_spectrum, resolution, \
                                        max_std_continuum = sigma, \
                                        continuum_model = sun_continuum_model, \
                                        max_continuum_diff=max_continuum_diff, \
                                        fixed_wave_step=fixed_wave_step)
    ispec.write_continuum_regions(sun_continuum_regions, "example_sun_fe_lines_continuum.txt")

def find_continuum_regions_in_segments():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Find continuum regions in segments ----------------------------------------
    logging.info("Finding continuum regions...")
    resolution = 80000
    sigma = 0.001
    max_continuum_diff = 1.0
    fixed_wave_step = 0.05
    # Limit the search to given segments
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    limited_sun_continuum_regions = ispec.find_continuum(sun_spectrum, resolution, \
                                            segments=segments, max_std_continuum = sigma, \
                                            continuum_model = sun_continuum_model, \
                                            max_continuum_diff=max_continuum_diff, \
                                            fixed_wave_step=fixed_wave_step)
    ispec.write_continuum_regions(limited_sun_continuum_regions, \
            "example_limited_sun_continuum_region.txt")


def find_linemasks():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Find linemasks ------------------------------------------------------------
    logging.info("Finding line masks...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

    resolution = 80000
    smoothed_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, resolution)
    min_depth = 0.05
    max_depth = 1.00
    vel_telluric = 17.79  # km/s
    sun_linemasks = ispec.find_linemasks(sun_spectrum, sun_continuum_model, \
                            atomic_linelist=atomic_linelist, \
                            max_atomic_wave_diff = 0.005, \
                            telluric_linelist=telluric_linelist, \
                            vel_telluric=vel_telluric, \
                            minimum_depth=min_depth, maximum_depth=max_depth, \
                            smoothed_spectrum=smoothed_sun_spectrum, \
                            check_derivatives=False, \
                            discard_gaussian=False, discard_voigt=True )
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (sun_linemasks['wave (nm)'] == 0)
    sun_linemasks = sun_linemasks[~rejected_by_atomic_line_not_found]

    ispec.write_line_regions(sun_linemasks, "example_sun_linemasks.txt")


def calculate_barycentric_velocity():
    mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")
    #--- Barycentric velocity correction from observation date/coordinates ---------
    logging.info("Calculating barycentric velocity correction...")
    day = 15
    month = 2
    year = 2012
    hours = 0
    minutes = 0
    seconds = 0
    ra_hours = 19
    ra_minutes = 50
    ra_seconds = 46.99
    dec_degrees = 8
    dec_minutes = 52
    dec_seconds = 5.96

    # Project velocity toward star
    barycentric_vel = ispec.calculate_barycentric_velocity_correction((year, month, day, \
                                    hours, minutes, seconds), (ra_hours, ra_minutes, \
                                    ra_seconds, dec_degrees, dec_minutes, dec_seconds))
    #--- Correcting barycentric velocity -------------------------------------------
    corrected_spectrum = ispec.correct_velocity(mu_cas_a_spectrum, barycentric_vel)
    return corrected_spectrum

def estimate_snr_from_flux():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Estimate SNR from flux ----------------------------------------------------
    logging.info("Estimating SNR...")
    num_points = 10
    estimated_snr = ispec.estimate_snr(sun_spectrum['flux'], num_points=num_points)
    return estimated_snr


def estimate_errors_from_snr():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Calculate errors based on SNR ---------------------------------------------
    snr = 100
    sun_spectrum['err'] = sun_spectrum['flux'] / snr
    return sun_spectrum


def clean_spectrum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Clean fluxes and errors ---------------------------------------------------
    logging.info("Cleaning fluxes and errors...")
    flux_base = 0.0
    flux_top = 1.0
    err_base = 0.0
    err_top = 1.0
    ffilter = (sun_spectrum['flux'] > flux_base) & (sun_spectrum['flux'] <= flux_top)
    efilter = (sun_spectrum['err'] > err_base) & (sun_spectrum['err'] <= err_top)
    wfilter = np.logical_and(ffilter, efilter)
    clean_sun_spectrum = sun_spectrum[wfilter]
    return clean_sun_spectrum


def clean_telluric_regions():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Clean regions that may be affected by tellurics ---------------------------
    logging.info("Cleaning tellurics...")

    telluric_lines_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"
    telluric_linelist = ispec.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)

    # - Filter regions that may be affected by telluric lines
    rv = 0.0
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(sun_spectrum['waveobs'], \
                                telluric_linelist[dfilter], min_velocity=-rv+min_vel, \
                                max_velocity=-rv+max_vel)
    clean_sun_spectrum = sun_spectrum[tfilter]
    return clean_sun_spectrum


def adjust_line_masks():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Adjust line masks ---------------------------------------------------------
    resolution = 80000
    smoothed_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, resolution)
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    linemasks = ispec.adjust_linemasks(smoothed_sun_spectrum, line_regions, max_margin=0.5)
    return linemasks

def create_segments_around_linemasks():
    #---Create segments around linemasks -------------------------------------------
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    segments = ispec.create_segments_around_lines(line_regions, margin=0.25)
    return segments


def determine_abundances_from_EW():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Fit lines -----------------------------------------------------------------
    logging.info("Fitting lines...")
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)


    vel_telluric = 17.79 # km/s
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.adjust_linemasks(sun_spectrum, line_regions, max_margin=0.5)
    # Spectrum should be already radial velocity corrected
    linemasks = ispec.fit_lines(line_regions, sun_spectrum, sun_continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                max_atomic_wave_diff = 0.005, \
                                telluric_linelist = telluric_linelist, \
                                check_derivatives = False, \
                                vel_telluric = vel_telluric, discard_gaussian=False, \
                                discard_voigt=True, free_mu=False)
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (linemasks['wave (nm)'] == 0)
    linemasks = linemasks[~rejected_by_atomic_line_not_found]
    # Exclude lines that may be affected by tellurics
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
    linemasks = linemasks[~rejected_by_telluric_line]

    #--- Determining abundances by EW of the previously fitted lines ---------------
    # Parameters
    teff = 5777.0
    logg = 4.44
    MH = 0.00
    microturbulence_vel = 1.0

    # Selected model amtosphere and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)
    spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(atmosphere_layers, \
            teff, logg, MH, linemasks, solar_abundances, microturbulence_vel = microturbulence_vel, \
            verbose=1)

    print "[X/H]: %.2f" % np.median(x_over_h)
    print "[X/Fe]: %.2f" % np.median(x_over_fe)

def synthesize_spectrum():
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    teff = 5777.0
    logg = 4.44
    MH = 0.00
    microturbulence_vel = 1.0
    macroturbulence = 0.0
    vsini = 2.0
    limb_darkening_coeff = 0.0
    resolution = 300000
    wave_step = 0.001

    # Wavelengths to synthesis
    #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    regions = None
    wave_base = 515.0 # Magnesium triplet region
    wave_top = 525.0

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    fixed_abundances = None # No fixed abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg


    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)

    # Synthesis
    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, linelist=atomic_linelist, abundances=solar_abundances, \
            fixed_abundances=fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1)
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    ispec.write_spectrum(synth_spectrum, "example_synth.s")
    return synth_spectrum


def add_noise_to_spectrum():
    """
    Add noise to an spectrum (ideally to a synthetic one) based on a given SNR.
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Adding poisson noise -----------------------------------------------------
    snr = 100
    distribution = "poisson" # "gaussian"
    noisy_sun_spectrum = ispec.add_noise(sun_spectrum, snr, distribution)
    return noisy_sun_spectrum


def determine_astrophysical_parameters_using_synth_spectra():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Modelize spectra ----------------------------------------------------------
    # Parameters
    initial_teff = 5750.0
    initial_logg = 4.5
    initial_MH = 0.00
    initial_vmic = 1.0
    initial_vmac = 0.0
    initial_vsini = 2.0
    initial_limb_darkening_coeff = 0.0
    initial_R = 80000
    max_iterations = 20

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    free_params = ["teff", "logg", "MH", "vmic", "vmac"]
    #free_params = ["R"]

    # Free individual element abundance
    free_abundances = None
    # - For letting free a given element abundance:
    #free_abundances = np.recarray((1, ), dtype=[('code', int),('Abund', float)])
    #element_abundance = 12 # Mg (Magnessium)
    #free_abundances['code'] = element_abundance
    #free_abundances['Abund'] = abundances['Abund'][abundances['code'] == int(element_abundance)]

    # Fe 1/2 regions
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")

    ### Add also regions from the wings of strong lines:
    # H beta
    hbeta_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Hbeta.txt")
    hbeta_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Hbeta_segments.txt")
    line_regions = np.hstack((line_regions, hbeta_lines))
    segments = np.hstack((segments, hbeta_segments))
    # H alpha
    halpha_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Halpha.txt")
    halpha_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Halpha_segments.txt")
    line_regions = np.hstack((line_regions, halpha_lines))
    segments = np.hstack((segments, halpha_segments))
    # Magnesium triplet
    mgtriplet_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_MgTriplet.txt")
    mgtriplet_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_MgTriplet_segments.txt")
    line_regions = np.hstack((line_regions, mgtriplet_lines))
    segments = np.hstack((segments, mgtriplet_segments))


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks = \
            ispec.modelize_spectrum(sun_spectrum, sun_continuum_model, \
            modeled_layers_pack, atomic_linelist, solar_abundances, free_abundances, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
            linemasks=line_regions, max_iterations=max_iterations)
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    ispec.write_spectrum(modeled_synth_spectrum, "example_modeled_synth.s")
    return obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks


def determine_astrophysical_parameters_from_EW():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Continuum fit -------------------------------------------------------------
    # One spline per each 5 nm
    model = "Splines" # "Polynomy"
    degree = 3
    nknots = None # Automatic: 1 spline every 1 nm
    from_resolution = 80000

    # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
    order='median+max'
    median_wave_range=0.01
    max_wave_range=1.0

    strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                                ignore=strong_lines, \
                                nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True)
    #--- Fit lines -----------------------------------------------------------------
    logging.info("Fitting lines...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)


    vel_telluric = 17.79 # km/s
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.adjust_linemasks(sun_spectrum, line_regions, max_margin=0.5)
    # Spectrum should be already radial velocity corrected
    linemasks = ispec.fit_lines(line_regions, sun_spectrum, sun_continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                max_atomic_wave_diff = 0.005, \
                                telluric_linelist = telluric_linelist, \
                                check_derivatives = False, \
                                vel_telluric = vel_telluric, discard_gaussian=False, \
                                discard_voigt=True, free_mu=False)
    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (linemasks['wave (nm)'] == 0)
    linemasks = linemasks[~rejected_by_atomic_line_not_found]
    # Exclude lines that may be affected by tellurics
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
    linemasks = linemasks[~rejected_by_telluric_line]

    #--- Modelize spectra from EW --------------------------------------------------
    # Parameters
    initial_teff = 5750.0
    initial_logg = 4.5
    initial_MH = 0.00
    initial_vmic = 1.0

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, initial_teff, initial_logg, initial_MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg

    # Reduced equivalent width
    # Filter too weak/strong lines
    # * Criteria presented in paper of GALA
    efilter = np.logical_and(linemasks['ewr'] >= -5.8, linemasks['ewr'] <= -4.65)
    # Filter high excitation potential lines
    # * Criteria from Eric J. Bubar "Equivalent Width Abundance Analysis In Moog"
    efilter = np.logical_and(efilter, linemasks['lower state (eV)'] <= 5.0)
    ## Filter also bad fits
    efilter = np.logical_and(efilter, linemasks['rms'] < 0.05)

    results = ispec.modelize_spectrum_from_EW(linemasks[efilter], modeled_layers_pack, atomic_linelist,\
                        solar_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, \
                        teff_elements=["Fe 1"], vmic_elements=["Fe 2"], \
                        adjust_model_metalicity=False, \
                        max_iterations=20)
    params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params = results
    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params


def calculate_theoretical_ew_and_depth():
    #--- Calculate theoretical equivalent widths and depths for a linelist ---------
    # Parameters
    teff = 5777.0
    logg = 4.44
    MH = 0.00
    microturbulence_vel = 1.0

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv3_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/475_685nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/475_685nm.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_noABO/845_895nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1_noABO/655_1020nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom/300_1100nm.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM/300_1000nm.lst"
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg


    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)

    # Synthesis
    #output_wave, output_code, output_ew, output_depth = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
    new_atomic_linelist = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
            teff, logg, MH, \
            atomic_linelist, solar_abundances, microturbulence_vel=microturbulence_vel, \
            verbose=1, gui_queue=None, timeout=900)
    ispec.write_atomic_linelist(new_atomic_linelist, linelist_filename="example_linelist.txt")
    return new_atomic_linelist


def add_noise_to_spectrum():
    """
    Add noise to an spectrum (ideally to a synthetic one) based on a given SNR.
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
    #--- Adding poisson noise -----------------------------------------------------
    snr = 100
    distribution = "poisson" # "gaussian"
    noisy_sun_spectrum = ispec.add_noise(sun_spectrum, snr, distribution)
    return noisy_sun_spectrum


if __name__ == '__main__':
    read_write_spectrum()
    convert_air_to_vacuum()
    #plot()
    cut_spectrum_from_range()
    cut_spectrum_from_segments()
    determine_radial_velocity_with_mask()
    determine_radial_velocity_with_template()
    correct_radial_velocity()
    determine_tellurics_shift()
    degrade_resolution()
    resample_spectrum()
    coadd_spectra()
    normalize_spectrum_using_continuum_regions()
    normalize_spectrum_in_segments()
    normalize_whole_spectrum_strategy2()
    normalize_whole_spectrum_strategy1()
    normalize_whole_spectrum_strategy1_ignoring_strong_lines()
    filter_cosmic_rays()
    find_continuum_regions()
    find_continuum_regions_in_segments()
    find_linemasks()
    calculate_barycentric_velocity()
    estimate_snr_from_flux()
    estimate_errors_from_snr()
    clean_spectrum()
    clean_telluric_regions()
    adjust_line_masks()
    create_segments_around_linemasks()
    determine_abundances_from_EW()
    synthesize_spectrum()
    add_noise_to_spectrum()
    determine_astrophysical_parameters_using_synth_spectra()
    determine_astrophysical_parameters_from_EW()
    calculate_theoretical_ew_and_depth()
    pass
