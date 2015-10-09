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
import multiprocessing
from multiprocessing import Pool

################################################################################
#--- iSpec directory -------------------------------------------------------------
ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
#ispec_dir = '/home/virtual/shared/iSpec/'
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
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    ispec.write_spectrum(sun_spectrum, "example_sun.s")
    return sun_spectrum



def convert_air_to_vacuum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Converting wavelengths from air to vacuum and viceversa -------------------
    sun_spectrum_vacuum = ispec.air_to_vacuum(sun_spectrum)
    sun_spectrum_air = ispec.vacuum_to_air(sun_spectrum_vacuum)
    return sun_spectrum_vacuum, sun_spectrum_air

def plot():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    #--- Plotting (requires graphical interface) -----------------------------------
    logging.info("Plotting...")
    ispec.plot_spectra([sun_spectrum, mu_cas_spectrum])
    ispec.show_histogram(sun_spectrum['flux'])

def cut_spectrum_from_range():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Cut -----------------------------------------------------------------------
    logging.info("Cutting...")

    # - Keep points between two given wavelengths
    wfilter = ispec.create_wavelength_filter(sun_spectrum, wave_base=480.0, wave_top=670.0)
    cutted_sun_spectrum = sun_spectrum[wfilter]
    return cutted_sun_spectrum

def cut_spectrum_from_segments():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Cut -----------------------------------------------------------------------
    logging.info("Cutting...")
    # Keep only points inside a list of segments
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    wfilter = ispec.create_wavelength_filter(sun_spectrum, regions=segments)
    cutted_sun_spectrum = sun_spectrum[wfilter]
    return cutted_sun_spectrum


def determine_radial_velocity_with_mask():
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    #--- Radial Velocity determination with linelist mask --------------------------
    logging.info("Radial velocity determination with linelist mask...")
    # - Read atomic data
    mask_file = ispec_dir + "input/linelists/CCF/Narval.Sun.370_1048nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Arcturus.372_926nm/mask.lst""
    #mask_file = ispec_dir + "input/linelists/CCF/Atlas.Sun.372_926nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.A0.350_1095nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.F0.360_698nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.G2.375_679nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K0.378_679nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.K5.378_680nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/HARPS_SOPHIE.M5.400_687nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/Synthetic.Sun.350_1100nm/mask.lst"
    #mask_file = ispec_dir + "input/linelists/CCF/VALD.Sun.300_1100nm/mask.lst"
    ccf_mask = ispec.read_linelist_mask(mask_file)

    models, ccf = ispec.cross_correlate_with_mask(mu_cas_spectrum, ccf_mask, \
                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                            velocity_step=1.0, mask_depth=0.01, \
                            fourier=False)

    # Number of models represent the number of components
    components = len(models)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    return rv, rv_err, components


def determine_radial_velocity_with_template():
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    #--- Radial Velocity determination with template -------------------------------
    logging.info("Radial velocity determination with template...")
    # - Read synthetic template
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Arcturus.372_926nm/template.txt.gz")
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Sun.372_926nm/template.txt.gz")
    template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.350_1100nm/template.txt.gz")

    models, ccf = ispec.cross_correlate_with_template(mu_cas_spectrum, template, \
                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                            velocity_step=1.0, fourier=False)

    # Number of models represent the number of components
    components = len(models)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    return rv, rv_err, components

def correct_radial_velocity():
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    #--- Radial Velocity correction ------------------------------------------------
    logging.info("Radial velocity correction...")
    rv = -96.40 # km/s
    mu_cas_spectrum = ispec.correct_velocity(mu_cas_spectrum, rv)


def determine_tellurics_shift_with_mask():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Telluric velocity shift determination from spectrum --------------------------
    logging.info("Telluric velocity shift determination...")
    # - Telluric
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    linelist_telluric = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    models, ccf = ispec.cross_correlate_with_mask(sun_spectrum, linelist_telluric, \
                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                            velocity_step=0.5, mask_depth=0.01, \
                            fourier = False,
                            only_one_peak = True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s
    return bv, bv_err

def determine_tellurics_shift_with_template():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Telluric velocity shift determination from spectrum --------------------------
    logging.info("Telluric velocity shift determination...")
    # - Read synthetic template
    template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Tellurics.350_1100nm/template.txt.gz")

    models, ccf = ispec.cross_correlate_with_template(sun_spectrum, template, \
                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                            velocity_step=0.5, fourier=False, \
                            only_one_peak = True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s
    return bv, bv_err



def degrade_resolution():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Resolution degradation ----------------------------------------------------
    logging.info("Resolution degradation...")
    from_resolution = 80000
    to_resolution = 40000
    convolved_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, to_resolution, \
                                                    from_resolution=from_resolution)
    return convolved_sun_spectrum

def smooth_spectrum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Smoothing spectrum (resolution will be affected) --------------------------
    logging.info("Smoothing spectrum...")
    resolution = 80000
    smoothed_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, resolution)
    return smoothed_sun_spectrum


def resample_spectrum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Resampling  --------------------------------------------------------------
    logging.info("Resampling...")
    wavelengths = np.arange(480.0, 670.0, 0.001)
    resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths, method="bessel", zero_edges=True)
    #resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths, method="linear", zero_edges=True)
    return resampled_sun_spectrum


def coadd_spectra():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    #--- Resampling and combining --------------------------------------------------
    logging.info("Resampling and comibining...")
    wavelengths = np.arange(480.0, 670.0, 0.001)
    resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths, zero_edges=True)
    resampled_mu_cas_spectrum = ispec.resample_spectrum(mu_cas_spectrum, wavelengths, zero_edges=True)
    # Coadd previously resampled spectra
    coadded_spectrum = ispec.create_spectrum_structure(resampled_sun_spectrum['waveobs'])
    coadded_spectrum['flux'] = resampled_sun_spectrum['flux'] + resampled_mu_cas_spectrum['flux']
    coadded_spectrum['err'] = np.sqrt(np.power(resampled_sun_spectrum['err'],2) + \
                                    np.power(resampled_mu_cas_spectrum['err'],2))
    return coadded_spectrum


def merge_spectra():
    #--- Mergin spectra ------------------------------------------------------------
    logging.info("Mergin spectra...")
    left_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    right_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
    merged_spectrum = np.hstack((left_spectrum, right_spectrum))
    return merged_spectrum


def normalize_spectrum_using_continuum_regions():
    """
    Consider only continuum regions for the fit, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                            automatic_strong_line_detection=True, \
                            strong_line_probability=0.5, \
                            use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_spectrum_in_segments():
    """
    Fit continuum in each segment independently, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                            automatic_strong_line_detection=True, \
                            strong_line_probability=0.5, \
                            use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy2():
    """
    Use the whole spectrum, strategy 'max+median'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy1():
    """
    Use the whole spectrum, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def normalize_whole_spectrum_strategy1_ignoring_prefixed_strong_lines():
    """
    Use the whole spectrum but ignoring some strong lines, strategy 'median+max'
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

    #--- Continuum normalization ---------------------------------------------------
    logging.info("Continuum normalization...")
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    return normalized_sun_spectrum, sun_continuum_model


def filter_cosmic_rays():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Filtering cosmic rays -----------------------------------------------------
    # Spectrum should be already normalized
    cosmics = ispec.create_filter_cosmic_rays(sun_spectrum, sun_continuum_model, \
                                            resampling_wave_step=0.001, window_size=15, \
                                            variation_limit=0.01)
    clean_sun_spectrum = sun_spectrum[~cosmics]
    return clean_sun_spectrum


def find_continuum_regions():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
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
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
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
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Find linemasks ------------------------------------------------------------
    logging.info("Finding line masks...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

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

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (sun_linemasks['ew'] == 0)
    sun_linemasks = sun_linemasks[~rejected_by_zero_ew]

    # Select only iron lines
    iron = sun_linemasks['element'] == "Fe 1"
    iron = np.logical_or(iron, sun_linemasks['element'] == "Fe 2")
    iron_sun_linemasks = sun_linemasks[iron]

    ispec.write_line_regions(sun_linemasks, "example_sun_linemasks.txt")
    ispec.write_line_regions(sun_linemasks, "example_sun_fe_linemasks.txt")


def calculate_barycentric_velocity():
    mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
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
    corrected_spectrum = ispec.correct_velocity(mu_cas_spectrum, barycentric_vel)
    return corrected_spectrum

def estimate_snr_from_flux():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    ## WARNING: To compare SNR estimation between different spectra, they should
    ##          be homogeneously sampled (consider a uniform re-sampling)
    #--- Estimate SNR from flux ----------------------------------------------------
    logging.info("Estimating SNR from fluxes...")
    num_points = 10
    estimated_snr = ispec.estimate_snr(sun_spectrum['flux'], num_points=num_points)
    return estimated_snr

def estimate_snr_from_err():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Estimate SNR from errors --------------------------------------------------
    logging.info("Estimating SNR from errors...")
    efilter = sun_spectrum['err'] > 0
    filtered_sun_spectrum = sun_spectrum[efilter]
    if len(filtered_sun_spectrum) > 1:
        estimated_snr = np.median(filtered_sun_spectrum['flux'] / filtered_sun_spectrum['err'])
    else:
        # All the errors are set to zero and we cannot calculate SNR using them
        estimated_snr = 0
    return estimated_snr


def estimate_errors_from_snr():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Calculate errors based on SNR ---------------------------------------------
    snr = 100
    sun_spectrum['err'] = sun_spectrum['flux'] / snr
    return sun_spectrum


def clean_spectrum():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
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
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Telluric velocity shift determination from spectrum --------------------------
    logging.info("Telluric velocity shift determination...")
    # - Telluric
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    linelist_telluric = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    models, ccf = ispec.cross_correlate_with_mask(sun_spectrum, linelist_telluric, \
                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                            velocity_step=0.5, mask_depth=0.01, \
                            fourier = False,
                            only_one_peak = True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s

    #--- Clean regions that may be affected by tellurics ---------------------------
    logging.info("Cleaning tellurics...")

    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    # - Filter regions that may be affected by telluric lines
    #bv = 0.0
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(sun_spectrum['waveobs'], \
                                telluric_linelist[dfilter], min_velocity=-bv+min_vel, \
                                max_velocity=-bv+max_vel)
    clean_sun_spectrum = sun_spectrum[~tfilter]
    return clean_sun_spectrum


def adjust_line_masks():
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
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

def fit_lines_and_determine_ew(use_ares=False):
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Fit lines -----------------------------------------------------------------
    logging.info("Fitting lines...")
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

    # Read
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)


    vel_telluric = 17.79 # km/s
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.adjust_linemasks(normalized_sun_spectrum, line_regions, max_margin=0.5)
    # Spectrum should be already radial velocity corrected
    linemasks = ispec.fit_lines(line_regions, normalized_sun_spectrum, sun_continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                #max_atomic_wave_diff = 0.005, \
                                max_atomic_wave_diff = 0.00, \
                                telluric_linelist = telluric_linelist, \
                                smoothed_spectrum = None, \
                                check_derivatives = False, \
                                vel_telluric = vel_telluric, discard_gaussian=False, \
                                discard_voigt=True, \
                                free_mu=True, crossmatch_with_mu=False, closest_match=True)
    # Discard lines that are not cross matched with the same original element stored in the note
    linemasks = linemasks[linemasks['element'] == line_regions['note']]

    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (linemasks['wave (nm)'] == 0)
    linemasks = linemasks[~rejected_by_atomic_line_not_found]

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (linemasks['ew'] == 0)
    linemasks = linemasks[~rejected_by_zero_ew]

    # Exclude lines that may be affected by tellurics
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
    linemasks = linemasks[~rejected_by_telluric_line]

    if use_ares:
        # Replace the measured equivalent widths by the ones computed by ARES
        old_linemasks = linemasks.copy()
        ### Different rejection parameters (check ARES papers):
        ##   - http://adsabs.harvard.edu/abs/2007A%26A...469..783S
        ##   - http://adsabs.harvard.edu/abs/2015A%26A...577A..67S
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="0.995", tmp_dir=None, verbose=0)
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="3;5764,5766,6047,6052,6068,6076", tmp_dir=None, verbose=0)
        snr = 50
        linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="%s" % (snr), tmp_dir=None, verbose=0)


    ew = linemasks['ew']
    ew_err = linemasks['ew_err']

    return linemasks, ew, ew_err


def synthesize_spectrum(code="spectrum"):
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
    #wave_base = 470
    ##wave_top = 570
    #wave_top = 680

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    isotopes = ispec.read_isotope_data(isotope_file)

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

    # Enhance alpha elements + CNO abundances following MARCS standard composition
    alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = ispec.determine_abundance_enchancements(MH)
    abundances = ispec.enhance_solar_abundances(solar_abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH, code=code)

    # Synthesis
    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, atomic_linelist, isotopes, abundances, \
            fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1,
            code=code, use_molecules=use_molecules)
    ##--- Save spectrum ------------------------------------------------------------
    logging.info("Saving spectrum...")
    if code == "turbospectrum":
        synth_filename = "example_synth_turbo.s"
    elif code == "moog":
        synth_filename = "example_synth_moog.s"
    else:
        synth_filename = "example_synth.s"
    ispec.write_spectrum(synth_spectrum, synth_filename)
    return synth_spectrum


def add_noise_to_spectrum():
    """
    Add noise to an spectrum (ideally to a synthetic one) based on a given SNR.
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Adding poisson noise -----------------------------------------------------
    snr = 100
    distribution = "poisson" # "gaussian"
    noisy_sun_spectrum = ispec.add_noise(sun_spectrum, snr, distribution)
    return noisy_sun_spectrum

def generate_new_random_realizations_from_spectrum():
    """
    Considering fluxes as mean values and errors as standard deviation, generate
    N new random spectra.
    """
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    number = 10
    spectra = ispec.random_realizations(sun_spectrum, number, distribution="poisson")
    return spectra

def precompute_synthetic_grid(code="spectrum"):
    if code == "turbospectrum":
        precomputed_grid_dir = "example_grid_turbo/"
    elif code == "moog":
        precomputed_grid_dir = "example_grid_moog/"
    else:
        precomputed_grid_dir = "example_grid/"

    # - Read grid ranges from file
    #ranges_filename = "input/grid/grid_ranges.txt"
    #ranges = ascii.read(ranges_filename, delimiter="\t")
    # - or define them directly here (example of only 2 reference points):
    ranges = np.recarray((2,),  dtype=[('Teff', int), ('logg', float), ('MH', float)])
    ranges['Teff'][0] = 5500
    ranges['logg'][0] = 4.5
    ranges['MH'][0] = 0.0
    ranges['Teff'][1] = 3500
    ranges['logg'][1] = 1.5
    ranges['MH'][1] = 0.0

    # Wavelengths
    initial_wave = 480.0
    final_wave = 670.0
    step_wave = 0.001
    wavelengths = np.arange(initial_wave, final_wave, step_wave)

    to_resolution = 80000 # Individual files will not be convolved but the grid will be (for fast comparison)
    number_of_processes = 1 # It can be parallelized for computers with multiple processors


    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    isotopes = ispec.read_isotope_data(isotope_file)

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


    ispec.precompute_synthetic_grid(precomputed_grid_dir, ranges, wavelengths, to_resolution, \
                                    modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, \
                                    segments=None, number_of_processes=number_of_processes, \
                                    code=code, use_molecules=use_molecules)


def determine_astrophysical_parameters_using_synth_spectra(code="spectrum"):
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    #normalized_sun_spectrum['flux'] *= 0.99
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Model spectra ----------------------------------------------------------
    # Parameters
    initial_teff = 5750.0
    initial_logg = 4.5
    initial_MH = 0.00
    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
    initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
    #initial_vsini = 2.0
    initial_vsini = 0.0
    initial_limb_darkening_coeff = 0.0
    initial_R = 80000
    max_iterations = 6

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    free_params = ["teff", "logg", "MH", "vmic", "vmac"]

    # Free individual element abundance
    free_abundances = None

    # Fe 1/2 regions
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_golden_lines.txt")
    line_regions = ispec.adjust_linemasks(normalized_sun_spectrum, line_regions, max_margin=0.5)
    # Read segments if we have them or...
    #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    # ... or we can create the segments on the fly:
    segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

    ### Add also regions from the wings of strong lines:
    ## H beta
    #hbeta_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Hbeta.txt")
    #hbeta_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Hbeta_segments.txt")
    #line_regions = np.hstack((line_regions, hbeta_lines))
    #segments = np.hstack((segments, hbeta_segments))
    ## H alpha
    #halpha_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Halpha.txt")
    #halpha_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Halpha_segments.txt")
    #line_regions = np.hstack((line_regions, halpha_lines))
    #segments = np.hstack((segments, halpha_segments))
    ## Magnesium triplet
    #mgtriplet_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_MgTriplet.txt")
    #mgtriplet_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_MgTriplet_segments.txt")
    #line_regions = np.hstack((line_regions, mgtriplet_lines))
    #segments = np.hstack((segments, mgtriplet_segments))

    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_sun_spectrum, sun_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
            linemasks=line_regions, \
            enhance_abundances=False, \
            use_errors = True, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code, use_molecules=use_molecules)
    ##--- Save results -------------------------------------------------------------
    logging.info("Saving results...")
    if code == "turbospectrum":
        dump_file = "example_results_synth_turbo.dump"
    elif code == "moog":
        dump_file = "example_results_synth_moog.dump"
    else:
        dump_file = "example_results_synth.dump"
    logging.info("Saving results...")
    ispec.save_results(dump_file, (params, errors, abundances_found, status, stats_linemasks))
    # If we need to restore the results from another script:
    params, errors, abundances_found, status, stats_linemasks = ispec.restore_results(dump_file)

    logging.info("Saving synthetic spectrum...")
    if code == "turbospectrum":
        synth_filename = "example_modeled_synth_turbo.s"
    elif code == "moog":
        synth_filename = "example_modeled_synth_moog.s"
    else:
        synth_filename = "example_modeled_synth.s"
    ispec.write_spectrum(modeled_synth_spectrum, synth_filename)

    return obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks


def determine_astrophysical_parameters_using_synth_spectra_and_precomputed_grid(code="spectrum"):
    ############################################################################
    # WARNING !!!
    #  This routine depends on the previous precomputation of the synthetic grid
    ############################################################################
    if code == "turbospectrum":
        precomputed_grid_dir = "example_grid_turbo/"
    elif code == "moog":
        precomputed_grid_dir = "example_grid_moog/"
    else:
        precomputed_grid_dir = "example_grid/"


    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")

    #--- Model spectra ----------------------------------------------------------
    # Parameters
    initial_R = 80000
    max_iterations = 6

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    free_params = ["teff", "logg", "MH", "vmic", "vmac"]

    # Free individual element abundance
    free_abundances = None

    # Fe 1/2 regions
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.adjust_linemasks(normalized_sun_spectrum, line_regions, max_margin=0.5)
    # Read segments if we have them or...
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    # ... or we can create the segments on the fly:
    #segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

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

    initial_teff, initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff = \
            ispec.estimate_initial_ap(normalized_sun_spectrum, precomputed_grid_dir, initial_R, line_regions)
    print "Initial estimation:", initial_teff, initial_logg, initial_MH, \
            initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff

    #--- Change LOG level ----------------------------------------------------------
    LOG_LEVEL = "warning"
    #LOG_LEVEL = "info"
    logger = logging.getLogger() # root logger, common for all
    logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_sun_spectrum, sun_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
            linemasks=line_regions, \
            enhance_abundances=True, \
            precomputed_grid_dir = precomputed_grid_dir, \
            use_errors = True, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code, use_molecules=use_molecules)

    #--- Change LOG level ----------------------------------------------------------
    #LOG_LEVEL = "warning"
    LOG_LEVEL = "info"
    logger = logging.getLogger() # root logger, common for all
    logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))

    ##--- Save results -------------------------------------------------------------
    logging.info("Saving results...")
    if code == "turbospectrum":
        dump_file = "example_results_synth_precomputed_turbo.dump"
    elif code == "moog":
        dump_file = "example_results_synth_precomputed_moog.dump"
    else:
        dump_file = "example_results_synth_precomputed.dump"
    logging.info("Saving results...")
    ispec.save_results(dump_file, (params, errors, abundances_found, status, stats_linemasks))
    # If we need to restore the results from another script:
    params, errors, abundances_found, status, stats_linemasks = ispec.restore_results(dump_file)

    logging.info("Saving synthetic spectrum...")
    if code == "turbospectrum":
        synth_filename = "example_modeled_synth_precomputed_turbo.s"
    elif code == "moog":
        synth_filename = "example_modeled_synth_precomputed_moog.s"
    else:
        synth_filename = "example_modeled_synth_precomputed.s"
    ispec.write_spectrum(modeled_synth_spectrum, synth_filename)
    return obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks




def determine_abundances_using_synth_spectra(code="spectrum"):
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Model spectra ----------------------------------------------------------
    # Parameters
    initial_teff = 5777.0
    initial_logg = 4.43
    initial_MH = 0.00
    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
    initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
    initial_vsini = 2.0
    initial_limb_darkening_coeff = 0.0
    initial_R = 80000
    max_iterations = 6

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)


    # Free parameters
    #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
    free_params = []

    # Free individual element abundance (WARNING: it should be coherent with the selecte line regions!)
    free_abundances = ispec.create_free_abundances_structure(["Fe"], chemical_elements, solar_abundances)

    # Fe 1/2 regions
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
    line_regions = ispec.adjust_linemasks(normalized_sun_spectrum, line_regions, max_margin=0.5)

    # Read segments if we have them or...
    segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
    # ... or we can create the segments on the fly:
    #segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

    obs_spec, modeled_synth_spectrum, params, errors, abundances_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_sun_spectrum, sun_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, initial_teff, \
            initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
            initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
            linemasks=line_regions, \
            enhance_abundances=True, \
            use_errors = True, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code, use_molecules=use_molecules)

    ##--- Save results -------------------------------------------------------------
    if code == "turbospectrum":
        dump_file = "example_results_synth_abundances_turbo.dump"
    elif code == "moog":
        dump_file = "example_results_synth_abundances_moog.dump"
    else:
        dump_file = "example_results_synth_abundances.dump"
    logging.info("Saving results...")
    ispec.save_results(dump_file, (params, errors, abundances_found, status, stats_linemasks))
    # If we need to restore the results from another script:
    params, errors, abundances_found, status, stats_linemasks = ispec.restore_results(dump_file)

    logging.info("Saving synthetic spectrum...")
    if code == "turbospectrum":
        synth_filename = "example_modeled_synth_abundances_turbo.s"
    elif code == "moog":
        synth_filename = "example_modeled_synth_abundances_moog.s"
    else:
        synth_filename = "example_modeled_synth_abundances.s"
    ispec.write_spectrum(modeled_synth_spectrum, synth_filename)

    return obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks


def determine_astrophysical_parameters_from_ew(code="spectrum", use_ares=False):
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Read lines and adjust them ------------------------------------------------
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines_biglist.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_golden_lines.txt")
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/GES/v5/atomic_lines_yy_fe.txt")
    line_regions = ispec.adjust_linemasks(sun_spectrum, line_regions, max_margin=0.5)
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Fit lines -----------------------------------------------------------------
    logging.info("Fitting lines...")

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

    # To reduce bad cross-matches, we reduce the linelist to those elements we are going to analyse
    efilter = None
    for element in np.unique(line_regions['note']):
        if efilter is None:
            efilter = atomic_linelist['element'] == element
        else:
            efilter = np.logical_or(efilter, atomic_linelist['element'] == element)
    atomic_linelist = atomic_linelist[efilter]


    vel_telluric = 17.79 # km/s
    #continuum_adjustment_margin = 0.05 # Allow +/-5% free baseline fit around continuum
    continuum_adjustment_margin = 0.0
    # Spectrum should be already radial velocity corrected
    linemasks = ispec.fit_lines(line_regions, normalized_sun_spectrum, sun_continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                max_atomic_wave_diff = 0.005, \
                                telluric_linelist = telluric_linelist, \
                                smoothed_spectrum = None, \
                                check_derivatives = False, \
                                vel_telluric = vel_telluric, discard_gaussian=False, \
                                discard_voigt=True, \
                                free_mu=True, crossmatch_with_mu=False, closest_match=True, \
                                continuum_adjustment_margin=continuum_adjustment_margin)
    # Discard lines that are not cross matched with the same original element stored in the note
    linemasks = linemasks[linemasks['element'] == line_regions['note']]

    # Discard bad masks
    flux_peak = normalized_sun_spectrum['flux'][linemasks['peak']]
    flux_base = normalized_sun_spectrum['flux'][linemasks['base']]
    flux_top = normalized_sun_spectrum['flux'][linemasks['top']]
    bad_mask = np.logical_or(linemasks['wave_peak'] <= linemasks['wave_base'], linemasks['wave_peak'] >= linemasks['wave_top'])
    bad_mask = np.logical_or(bad_mask, flux_peak >= flux_base)
    bad_mask = np.logical_or(bad_mask, flux_peak >= flux_top)
    linemasks = linemasks[~bad_mask]

    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (linemasks['wave (nm)'] == 0)
    linemasks = linemasks[~rejected_by_atomic_line_not_found]

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (linemasks['ew'] == 0)
    linemasks = linemasks[~rejected_by_zero_ew]

    # Exclude lines that may be affected by tellurics
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
    linemasks = linemasks[~rejected_by_telluric_line]

    if use_ares:
        # Replace the measured equivalent widths by the ones computed by ARES
        old_linemasks = linemasks.copy()
        ### Different rejection parameters (check ARES papers):
        ##   - http://adsabs.harvard.edu/abs/2007A%26A...469..783S
        ##   - http://adsabs.harvard.edu/abs/2015A%26A...577A..67S
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="0.995", tmp_dir=None, verbose=0)
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="3;5764,5766,6047,6052,6068,6076", tmp_dir=None, verbose=0)
        snr = 50
        linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="%s" % (snr), tmp_dir=None, verbose=0)


    #--- Model spectra from EW --------------------------------------------------
    # Parameters
    initial_teff = 5750.0
    initial_logg = 4.5
    initial_MH = 0.00
    initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
    max_iterations = 10

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"


    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

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
    #efilter = np.logical_and(linemasks['ewr'] >= -5.8, linemasks['ewr'] <= -4.65)
    efilter = np.logical_and(linemasks['ewr'] >= -6.0, linemasks['ewr'] <= -4.3)
    # Filter high excitation potential lines
    # * Criteria from Eric J. Bubar "Equivalent Width Abundance Analysis In Moog"
    efilter = np.logical_and(efilter, linemasks['lower state (eV)'] <= 5.0)
    efilter = np.logical_and(efilter, linemasks['lower state (eV)'] >= 0.5)
    ## Filter also bad fits
    efilter = np.logical_and(efilter, linemasks['rms'] < 0.05)

    if code in ["moog", "width"]:
        adjust_model_metalicity = True
    else:
        adjust_model_metalicity = False

    results = ispec.model_spectrum_from_ew(linemasks[efilter], modeled_layers_pack, atomic_linelist,\
                        solar_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, \
                        free_params=["teff", "logg", "vmic"], \
                        adjust_model_metalicity=adjust_model_metalicity, \
                        max_iterations=max_iterations, \
                        enhance_abundances=True, \
                        outliers_detection = "robust", \
                        outliers_weight_limit = 0.90, \
                        tmp_dir = None, \
                        code=code)
    params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params = results

    ##--- Save results -------------------------------------------------------------
    logging.info("Saving results...")
    if code == "turbospectrum":
        dump_file = "example_results_ew_turbo.dump"
    elif code == "moog":
        dump_file = "example_results_ew_moog.dump"
    else:
        dump_file = "example_results_ew.dump"

    ispec.save_results(dump_file, (params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params))
    # If we need to restore the results from another script:
    params, errors, status, x_over_h, selected_x_over_h, fitted_lines_param = ispec.restore_results(dump_file)

    return params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params


def determine_abundances_from_ew(code="spectrum", use_ares=False):
    sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
    #--- Read lines and adjust them ------------------------------------------------
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines_biglist.txt")
    #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_golden_lines.txt")
    line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/GES/v5/atomic_lines_yy_fe.txt")
    line_regions = ispec.adjust_linemasks(sun_spectrum, line_regions, max_margin=0.5)
    #--- Continuum fit -------------------------------------------------------------
    model = "Splines" # "Polynomy"
    degree = 2
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
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)
    #--- Normalize -------------------------------------------------------------
    normalized_sun_spectrum = ispec.normalize_spectrum(sun_spectrum, sun_continuum_model, consider_continuum_errors=False)
    # Use a fixed value because the spectrum is already normalized
    sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")
    #--- Fit lines -----------------------------------------------------------------
    logging.info("Fitting lines...")

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

    # Load chemical information and linelist
    if code == "turbospectrum":
        use_molecules = True # Only for turbo
    else:
        use_molecules = False
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    if code == "turbospectrum":
        atomic_linelist_file = ispec_dir + "input/linelists/turbospectrum/GESv5_atom_hfs_iso.420_920nm/atomic_lines.bsyn"
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, code=code) ## TurboSpectrum
    else:
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements=chemical_elements, molecules=molecules, code=code)
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

    # To reduce bad cross-matches, we reduce the linelist to those elements we are going to analyse
    efilter = None
    for element in np.unique(line_regions['note']):
        if efilter is None:
            efilter = atomic_linelist['element'] == element
        else:
            efilter = np.logical_or(efilter, atomic_linelist['element'] == element)
    atomic_linelist = atomic_linelist[efilter]


    vel_telluric = 17.79 # km/s
    # Spectrum should be already radial velocity corrected
    linemasks = ispec.fit_lines(line_regions, normalized_sun_spectrum, sun_continuum_model, \
                                atomic_linelist = atomic_linelist, \
                                max_atomic_wave_diff = 0.005, \
                                telluric_linelist = telluric_linelist, \
                                smoothed_spectrum = None, \
                                check_derivatives = False, \
                                vel_telluric = vel_telluric, discard_gaussian=False, \
                                discard_voigt=True, \
                                free_mu=True, crossmatch_with_mu=False, closest_match=True)
    # Discard lines that are not cross matched with the same original element stored in the note
    linemasks = linemasks[linemasks['element'] == line_regions['note']]

    # Exclude lines that have not been successfully cross matched with the atomic data
    # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
    rejected_by_atomic_line_not_found = (linemasks['wave (nm)'] == 0)
    linemasks = linemasks[~rejected_by_atomic_line_not_found]

    # Exclude lines with EW equal to zero
    rejected_by_zero_ew = (linemasks['ew'] == 0)
    linemasks = linemasks[~rejected_by_zero_ew]

    # Exclude lines that may be affected by tellurics
    rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)
    linemasks = linemasks[~rejected_by_telluric_line]

    if use_ares:
        # Replace the measured equivalent widths by the ones computed by ARES
        old_linemasks = linemasks.copy()
        ### Different rejection parameters (check ARES papers):
        ##   - http://adsabs.harvard.edu/abs/2007A%26A...469..783S
        ##   - http://adsabs.harvard.edu/abs/2015A%26A...577A..67S
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="0.995", tmp_dir=None, verbose=0)
        #linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="3;5764,5766,6047,6052,6068,6076", tmp_dir=None, verbose=0)
        snr = 50
        linemasks = ispec.update_ew_with_ares(normalized_sun_spectrum, linemasks, rejt="%s" % (snr), tmp_dir=None, verbose=0)


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
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH, code=code)
    spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(atmosphere_layers, \
            teff, logg, MH, linemasks, solar_abundances, microturbulence_vel = microturbulence_vel, \
            verbose=1, code=code)

    bad = np.isnan(x_over_h)
    fe1 = linemasks['element'] == "Fe 1"
    fe2 = linemasks['element'] == "Fe 2"
    print "[Fe 1/H]: %.2f" % np.median(x_over_h[np.logical_and(fe1, ~bad)])
    #print "[X/Fe]: %.2f" % np.median(x_over_fe[np.logical_and(fe1, ~bad)])
    print "[Fe 2/H]: %.2f" % np.median(x_over_h[np.logical_and(fe2, ~bad)])
    #print "[X/Fe]: %.2f" % np.median(x_over_fe[np.logical_and(fe2, ~bad)])



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

    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/VALD_atom.300_1100nm/atomic_lines.lst"
    atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_hfs_iso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_hfs_iso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/GESv4_atom_nohfs_noiso.475_685nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SEPv1.655_1020nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/Kurucz_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/NIST_atom.300_1100nm/atomic_lines.lst"
    #atomic_linelist_file = ispec_dir + "/input/linelists/SPECTRUM/SPECTRUM.300_1000nm/atomic_lines.lst"

    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
    molecules_file = ispec_dir + "/input/abundances/molecular_symbols.dat"

    solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    molecules = ispec.read_molecular_symbols(molecules_file)
    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, chemical_elements, molecules)
    isotopes = ispec.read_isotope_data(isotope_file)

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg

    # Enhance alpha elements + CNO abundances following MARCS standard composition
    alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = ispec.determine_abundance_enchancements(MH)
    abundances = ispec.enhance_solar_abundances(solar_abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)

    # Synthesis
    #output_wave, output_code, output_ew, output_depth = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
    new_atomic_linelist = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
            teff, logg, MH, \
            atomic_linelist, isotopes, abundances, microturbulence_vel=microturbulence_vel, \
            verbose=1, gui_queue=None, timeout=900)
    ispec.write_atomic_linelist(new_atomic_linelist, linelist_filename="example_linelist.txt")
    return new_atomic_linelist



def analyze(text, number):
    #--- Fake function to be used in the parallelization example -------------------
    multiprocessing.current_process().daemon=False
    import time

    # Print some text and wait for some seconds to finish
    print "Starting", text
    time.sleep(2 + number)
    print "... end of", number

def paralelize_code():
    number_of_processes = 2
    pool = Pool(number_of_processes)
    #--- Send 5 analyze processes to the pool which will execute 2 in parallel -----
    pool.apply_async(analyze, ["one", 1])
    pool.apply_async(analyze, ["two", 2])
    pool.apply_async(analyze, ["three", 3])
    pool.apply_async(analyze, ["four", 4])
    pool.apply_async(analyze, ["five", 5])
    pool.close()
    pool.join()


def estimate_vmic_from_empirical_relation():
    teff = 5500
    logg = 4.5
    MH = 0.0
    vmic = ispec.estimate_vmic(teff, logg, MH)
    print "VMIC:", vmic

def estimate_vmac_from_empirical_relation():
    teff = 5500
    logg = 4.5
    MH = 0.0
    vmac = ispec.estimate_vmac(teff, logg, MH)
    print "VMAC:", vmac

def generate_and_plot_YY_isochrone():
    import isochrones
    import matplotlib.pyplot as plt

    logage = 9.409
    age = np.power(10, logage) / 1e9 # Gyrs
    MH = 0.0 # [M/H] (dex)
    isochrone = isochrones.interpolate_isochrone(ispec_dir, age, MH)

    plt.plot(np.power(10, isochrone['logT']),  isochrone['logg'], marker='', ls='-', color="blue", label="[M/H] %.2f, %.2f Gyrs" % (MH, age))
    plt.xlabel("$T_{eff}$ (K)")
    plt.ylabel("$log(g)$ (dex)")
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend()
    plt.show()



def interpolate_atmosphere(code="spectrum"):
    #--- Synthesizing spectrum -----------------------------------------------------
    # Parameters
    teff = 5777.0
    logg = 4.44
    MH = 0.05

    # Selected model amtosphere, linelist and solar abundances
    #model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
    model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
    #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Validate parameters
    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of theatmospheric models."
        print msg

    # Prepare atmosphere model
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH, code=code)
    atmosphere_layers_file = "example_atmosphere.txt"
    atmosphere_layers_file = ispec.write_atmosphere(atmosphere_layers, teff, logg, MH, atmosphere_filename=atmosphere_layers_file, code=code)


if __name__ == '__main__':
    read_write_spectrum()
    convert_air_to_vacuum()
    #plot()
    cut_spectrum_from_range()
    cut_spectrum_from_segments()
    determine_radial_velocity_with_mask()
    determine_radial_velocity_with_template()
    correct_radial_velocity()
    determine_tellurics_shift_with_mask()
    determine_tellurics_shift_with_template()
    degrade_resolution()
    smooth_spectrum()
    resample_spectrum()
    coadd_spectra()
    merge_spectra()
    normalize_spectrum_using_continuum_regions()
    normalize_spectrum_in_segments()
    normalize_whole_spectrum_strategy2()
    normalize_whole_spectrum_strategy1()
    normalize_whole_spectrum_strategy1_ignoring_prefixed_strong_lines()
    filter_cosmic_rays()
    find_continuum_regions()
    find_continuum_regions_in_segments()
    find_linemasks()
    fit_lines_and_determine_ew(use_ares=False)
    fit_lines_and_determine_ew(use_ares=True)
    calculate_barycentric_velocity()
    estimate_snr_from_flux()
    estimate_snr_from_err()
    estimate_errors_from_snr()
    clean_spectrum()
    clean_telluric_regions()
    adjust_line_masks()
    create_segments_around_linemasks()
    synthesize_spectrum(code="spectrum")
    synthesize_spectrum(code="turbospectrum")
    synthesize_spectrum(code="moog")
    add_noise_to_spectrum()
    generate_new_random_realizations_from_spectrum()
    #precompute_synthetic_grid(code="spectrum")
    #precompute_synthetic_grid(code="turbospectrum")
    determine_astrophysical_parameters_using_synth_spectra(code="spectrum")
    determine_astrophysical_parameters_using_synth_spectra(code="turbospectrum")
    determine_astrophysical_parameters_using_synth_spectra(code="moog")
    #determine_astrophysical_parameters_using_synth_spectra_and_precomputed_grid(code="spectrum")
    #determine_astrophysical_parameters_using_synth_spectra_and_precomputed_grid(code="turbospectrum")
    #determine_astrophysical_parameters_using_synth_spectra_and_precomputed_grid(code="moog")
    determine_abundances_using_synth_spectra(code="spectrum")
    determine_abundances_using_synth_spectra(code="turbospectrum")
    determine_abundances_using_synth_spectra(code="moog")
    #determine_astrophysical_parameters_from_ew(code="spectrum", use_ares=False)
    #determine_astrophysical_parameters_from_ew(code="turbospectrum", use_ares=False)
    determine_astrophysical_parameters_from_ew(code="moog", use_ares=False)
    determine_astrophysical_parameters_from_ew(code="width", use_ares=False)
    #determine_astrophysical_parameters_from_ew(code="spectrum", use_ares=True)
    #determine_astrophysical_parameters_from_ew(code="turbospectrum", use_ares=True)
    determine_astrophysical_parameters_from_ew(code="moog", use_ares=True)
    determine_astrophysical_parameters_from_ew(code="width", use_ares=True)
    #determine_abundances_from_ew(code="spectrum", use_ares=False)
    #determine_abundances_from_ew(code="turbospectrum", use_ares=False)
    determine_abundances_from_ew(code="moog", use_ares=False)
    determine_abundances_from_ew(code="width", use_ares=False)
    #determine_abundances_from_ew(code="spectrum", use_ares=True)
    #determine_abundances_from_ew(code="turbospectrum", use_ares=True)
    determine_abundances_from_ew(code="moog", use_ares=True)
    determine_abundances_from_ew(code="width", use_ares=True)
    calculate_theoretical_ew_and_depth()
    paralelize_code()
    estimate_vmic_from_empirical_relation()
    estimate_vmac_from_empirical_relation()
    #generate_and_plot_YY_isochrone()
    interpolate_atmosphere(code="spectrum")
    interpolate_atmosphere(code="turbospectrum")
    interpolate_atmosphere(code="moog")
    pass
