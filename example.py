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
from astropy.io import ascii

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


#--- Reading spectra -----------------------------------------------------------
logging.info("Reading spectra")
sun_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_sun.s.gz")
mu_cas_a_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/narval_mu_cas.s.gz")


#--- Converting wavelengths from air to vacuum and viceversa -------------------
sun_spectrum_vacuum = ispec.air_to_vacuum(sun_spectrum)
sun_spectrum_air = ispec.vacuum_to_air(sun_spectrum_vacuum)


#--- Plotting (requires graphical interface) -----------------------------------
logging.info("Plotting...")
ispec.plot_spectra([sun_spectrum, mu_cas_a_spectrum])
ispec.show_histogram(sun_spectrum['flux'])


#--- Cut -----------------------------------------------------------------------
logging.info("Cutting...")

# - Keep points between two given wavelengths
wfilter = ispec.create_wavelength_filter(sun_spectrum, wave_base=480.0, wave_top=670.0)
cutted_sun_spectrum = sun_spectrum[wfilter]

# - Keep only points inside a list of segments
segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
wfilter = ispec.create_wavelength_filter(sun_spectrum, regions=segments)
cutted_sun_spectrum = sun_spectrum[wfilter]



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
ccf_mask = ascii.read(mask_file)._data

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


#--- Radial Velocity correction ------------------------------------------------
logging.info("Radial velocity correction...")
mu_cas_a_spectrum = ispec.correct_velocity(mu_cas_a_spectrum, rv)


#--- Barycentric Velocity determination from spectrum --------------------------
logging.info("Barycentric velocity determination...")
# - Telluric
telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"
linelist_telluric = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

xcoord, fluxes, errors = ispec.build_velocity_profile(sun_spectrum, \
                                        linelist=linelist_telluric, lower_velocity_limit=-100.0, \
                                        upper_velocity_limit=100.0, velocity_step=0.5)

models = ispec.modelize_velocity_profile(xcoord, fluxes, errors, only_one_peak=True)
bv = np.round(models[0].mu(), 2) # km/s


#--- Resolution degradation ----------------------------------------------------
logging.info("Resolution degradation...")
from_resolution = 80000
to_resolution = 40000
convolved_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, to_resolution, \
                                                from_resolution=from_resolution)
convolved_mu_cas_a_spectrum = ispec.convolve_spectrum(mu_cas_a_spectrum, to_resolution, \
                                                from_resolution=from_resolution)


#--- Resampling and combining --------------------------------------------------
logging.info("Resampling and comibining...")
wavelengths = np.arange(480.0, 670.0, 0.001)
resampled_sun_spectrum = ispec.resample_spectrum(sun_spectrum, wavelengths)
resampled_mu_cas_a_spectrum = ispec.resample_spectrum(mu_cas_a_spectrum, wavelengths)

# Coadd previously resampled spectra
total_wavelengths = len(resampled_sun_spectrum)
coadded_spectrum = ispec.create_spectrum_structure(resampled_sun_spectrum['waveobs'])
coadded_spectrum['flux'] = resampled_sun_spectrum['flux'] + resampled_mu_cas_a_spectrum['flux']
coadded_spectrum['err'] = np.sqrt(np.power(resampled_sun_spectrum['err'],2) + \
                                np.power(resampled_mu_cas_a_spectrum['err'],2))


#--- Fit continuum -------------------------------------------------------------
logging.info("Fit continuum...")

# EXAMPLE 1: Use a fixed value (useful when the spectrum is already normalized)
# =========
sun_continuum_model = ispec.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")


# EXAMPLE 2: Consider only continuum regions for the fit, strategy 'median+max'
# =========

# One spline per each 5 nm
model = "Splines" # "Polynomy"
degree = 3
nknots = None # Automatic: 1 spline every 5 nm
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

# EXAMPLE 3: Fit continuum in each segment independently, strategy 'median+max'
#==========
# One spline per each 5 nm
model = "Splines" # "Polynomy"
degree = 3
nknots = None # Automatic: 1 spline every 5 nm
from_resolution = 80000

# Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
order='median+max'
median_wave_range=0.01
max_wave_range=1.0

segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
sun_continuum_model = ispec.fit_continuum(sun_spectrum, from_resolution=from_resolution, \
                        independent_regions=segments, nknots=1, degree=degree,\
                        median_wave_range=median_wave_range, \
                        max_wave_range=max_wave_range, \
                        model=model, order=order, \
                        automatic_strong_line_detection=True)



# EXAMPLE 4: Use the whole spectrum, strategy 'max+median'
#==========
# One spline per each 5 nm
model = "Splines" # "Polynomy"
degree = 3
nknots = None # Automatic: 1 spline every 5 nm
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


# EXAMPLE 5: Use the whole spectrum, strategy 'median+max'
#==========
# One spline per each 5 nm
model = "Splines" # "Polynomy"
degree = 3
nknots = None # Automatic: 1 spline every 5 nm
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


# EXAMPLE 6: Use the whole spectrum but ignoring some strong lines, strategy 'median+max'
#==========
# One spline per each 5 nm
model = "Splines" # "Polynomy"
degree = 3
nknots = None # Automatic: 1 spline every 5 nm
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
normalized_sun_spectrum = ispec.create_spectrum_structure(sun_spectrum['waveobs'])
normalized_sun_spectrum['flux'] = sun_spectrum['flux'] \
                                    / sun_continuum_model(sun_spectrum['waveobs'])
normalized_sun_spectrum['err'] = sun_spectrum['err'] \
                                    / sun_continuum_model(sun_spectrum['waveobs'])


#--- Filtering cosmic rays -----------------------------------------------------
# Spectrum should be already normalized
cosmics = ispec.create_filter_cosmic_rays(sun_spectrum, sun_continuum_model, \
                                        resampling_wave_step=0.001, window_size=15, \
                                        variation_limit=0.01)
clean_sun_spectrum = sun_spectrum[~cosmics]


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

# Or limit the search to given segments
segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
limited_sun_continuum_regions = ispec.find_continuum(sun_spectrum, resolution, \
                                        segments=segments, max_std_continuum = sigma, \
                                        continuum_model = sun_continuum_model, \
                                        max_continuum_diff=max_continuum_diff, \
                                        fixed_wave_step=fixed_wave_step)
ispec.write_continuum_regions(limited_sun_continuum_regions, \
        "example_limited_sun_continuum_region.txt")


#--- Find linemasks ------------------------------------------------------------
logging.info("Finding line masks...")
sun_continuum_model = ispec.fit_continuum(sun_spectrum)
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
resolution = 80000
smoothed_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, resolution)
min_depth = 0.05
max_depth = 1.00
vel_telluric = 17.79  # km/s
sun_linemasks = ispec.find_linemasks(sun_spectrum, sun_continuum_model, \
                        atomic_linelist_file=atomic_linelist_file, \
                        chemical_elements_file=chemical_elements_file, \
                        molecules_file=molecules_file, \
                        max_atomic_wave_diff = 0.0005, \
                        telluric_linelist_file=telluric_linelist_file, \
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
corrected_spectrum = ispec.correct_velocity(mu_cas_a_spectrum, barycentric_vel)


#--- Estimate SNR from flux ----------------------------------------------------
logging.info("Estimating SNR...")
num_points = 10
estimated_snr = ispec.estimate_snr(sun_spectrum['flux'], num_points=num_points)


#--- Calculate errors based on SNR ---------------------------------------------
snr = 100
sun_spectrum['err'] = sun_spectrum['flux'] / snr


#--- Clean fluxes and errors ---------------------------------------------------
logging.info("Cleaning fluxes and errors...")
flux_base = 0.0
flux_top = 1.0
err_base = 0.0
err_top = 1.0
ffilter = (sun_spectrum['flux'] > flux_base) & (sun_spectrum['flux'] <= flux_top)
efilter = (sun_spectrum['err'] > err_base) & (sun_spectrum['err'] <= err_top)
wfilter = np.logical_and(ffilter, efilter)
sun_spectrum = sun_spectrum[wfilter]


#--- Clean regions that may be affected by tellurics ---------------------------
logging.info("Cleaning tellurics...")

telluric_lines_file = ispec_dir + "/input/linelists/CCF/Tellurics.standard.atm_air_model.txt"
linelist_telluric = ispec.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)

# - Filter regions that may be affected by telluric lines
rv = 0.0
min_vel = -30.0
max_vel = +30.0
# Only the 25% of the deepest ones:
dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
tfilter = ispec.create_filter_for_regions_affected_by_tellurics(sun_spectrum['waveobs'], \
                            linelist_telluric[dfilter], min_velocity=-rv+min_vel, \
                            max_velocity=-rv+max_vel)
clean_sun_spectrum = sun_spectrum[tfilter]


#--- Adjust line masks ---------------------------------------------------------
resolution = 80000
smoothed_sun_spectrum = ispec.convolve_spectrum(sun_spectrum, resolution)
line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
linemasks = ispec.adjust_linemasks(smoothed_sun_spectrum, line_regions, max_margin=0.5)
segments = ispec.create_segments_around_lines(linemasks, margin=0.25)


#---Create segments around linemasks -------------------------------------------
line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
segments = ispec.create_segments_around_lines(line_regions, margin=0.25)


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
vel_telluric = 17.79 # km/s
line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
line_regions = ispec.adjust_linemasks(sun_spectrum, line_regions, max_margin=0.5)
# Spectrum should be already radial velocity corrected
linemasks = ispec.fit_lines(line_regions, sun_spectrum, sun_continuum_model, \
                            atomic_linelist_file = atomic_linelist_file, \
                            chemical_elements_file = chemical_elements_file, \
                            molecules_file = molecules_file, \
                            max_atomic_wave_diff = 0.005, \
                            telluric_linelist_file = telluric_linelist_file, \
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

# Selected model amtosphere, linelist and solar abundances
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
solar_abundances = ispec.read_SPECTRUM_abundances(solar_abundances_file)

# Validate parameters
if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
    msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
            fall out of theatmospheric models."
    print msg

# Prepare atmosphere model
atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)
spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(atmosphere_layers, \
        teff, logg, MH, linemasks, solar_abundances, microturbulence_vel = 2.0, verbose=1)

print "[X/H]: %.2f" % np.median(x_over_h)
print "[X/Fe]: %.2f" % np.median(x_over_fe)


#--- Synthesizing spectra ------------------------------------------------------
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

solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

# Load model atmospheres
modeled_layers_pack = ispec.load_modeled_layers_pack(model)
# Load SPECTRUM linelist
linelist = ispec.read_SPECTRUM_linelist(atomic_linelist_file)
# Load SPECTRUM abundances
fixed_abundances = None # No fixed abundances
solar_abundances = ispec.read_SPECTRUM_abundances(solar_abundances_file)

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
        atmosphere_layers, teff, logg, MH, linelist=linelist, abundances=solar_abundances, \
        fixed_abundances=fixed_abundances, microturbulence_vel = microturbulence_vel, \
        macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
        R=resolution, regions=regions, verbose=1)


#--- Adding gaussian noise -----------------------------------------------------
snr = 100
distribution = "poisson" # "gaussian"
synth_spectrum = ispec.add_noise(synth_spectrum, snr, distribution)


#--- Modelize spectra ----------------------------------------------------------
# Parameters
initial_teff = 5777.0
initial_logg = 4.44
initial_MH = 0.00
initial_vmic = 1.0
initial_vmac = 0.0
initial_vsini = 2.0
initial_limb_darkening_coeff = 0.0
initial_R = 300000
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

solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
#solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

# Free parameters
#free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "limb_darkening_coeff"]
free_params = ["R"]

# Free individual element abundance
free_abundances = None
# - For letting free a given element abundance:
#free_abundances = np.recarray((1, ), dtype=[('code', int),('Abund', float)])
#element_abundance = 12 # Mg (Magnessium)
#free_abundances['code'] = element_abundance
#free_abundances['Abund'] = abundances['Abund'][abundances['code'] == int(element_abundance)]

# Regions
segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")

# Load model atmospheres
modeled_layers_pack = ispec.load_modeled_layers_pack(model)
# Load SPECTRUM linelist
linelist = ispec.read_SPECTRUM_linelist(atomic_linelist_file)
# Load SPECTRUM abundances
fixed_abundances = None # No fixed abundances
solar_abundances = ispec.read_SPECTRUM_abundances(solar_abundances_file)

obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks = \
        ispec.modelize_spectrum(sun_spectrum, sun_continuum_model, \
        modeled_layers_pack, linelist, solar_abundances, free_abundances, initial_teff, \
        initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
        initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
        linemasks=line_regions, max_iterations=max_iterations)

#--- Modelize spectra from EW --------------------------------------------------
# Parameters
initial_teff = 5750.0
initial_logg = 4.5
initial_MH = 0.00
initial_vmic = 0.8

# Validate parameters
if not ispec.valid_atmosphere_target(modeled_layers_pack, initial_teff, initial_logg, initial_MH):
    msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
            fall out of theatmospheric models."
    print msg

# Reduced equivalent width
ewr = np.log10(linemasks['ew']/linemasks['wave_peak'])
# Filter too weak/strong lines
# Filter criteria exposed in paper of GALA
efilter = np.logical_and(ewr >= -5.8, ewr <= -4.65)

results = ispec.modelize_spectrum_from_EW(linemasks[efilter], modeled_layers_pack, linelist,\
                    solar_abundances, initial_teff, initial_logg, initial_MH, initial_vmic, \
                    teff_elements=["Fe 1"], vmic_elements=["Fe 2"], max_iterations=20)
params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params = results


##--- Save spectrum ------------------------------------------------------------
logging.info("Saving spectrum...")
ispec.write_spectrum(sun_spectrum, "example_sun.s")
ispec.write_spectrum(mu_cas_a_spectrum, "example_mu_cas_a.s")
ispec.write_spectrum(synth_spectrum, "example_synth.s")
ispec.write_spectrum(modeled_synth_spectrum, "example_modeled_synth.s")


