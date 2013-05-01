#!/usr/bin/env python
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
import numpy as np
import logging

################################################################################
#--- SVE directory -------------------------------------------------------------
#sve_dir = '/home/blanco/shared/sve/'
sve_dir = './'
sys.path.insert(0, os.path.abspath(sve_dir))
import sve


#--- Change LOG level ----------------------------------------------------------
#LOG_LEVEL = "warning"
LOG_LEVEL = "info"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))
################################################################################


#--- Reading spectra -----------------------------------------------------------
logging.info("Reading spectra")
sun_spectrum = sve.read_spectrum(sve_dir + "/input/spectra/examples/narval_sun.s.gz")
mu_cas_a_spectrum = sve.read_spectrum(sve_dir + "/input/spectra/examples/narval_mu_cas.s.gz")


#--- Plotting (requires graphical interface) -----------------------------------
logging.info("Plotting...")
sve.plot_spectra([sun_spectrum, mu_cas_a_spectrum])
sve.show_histogram(sun_spectrum['flux'])


#--- Cut -----------------------------------------------------------------------
logging.info("Cutting...")

# - Keep points between two given wavelengths
wfilter = sve.create_wavelength_filter(sun_spectrum, wave_base=480.0, wave_top=670.0)
cutted_sun_spectrum = sun_spectrum[wfilter]

# - Keep only points inside a list of segments
segments = sve.read_segment_regions(sve_dir + "/input/regions/fe_lines_segments.txt")
wfilter = sve.create_wavelength_filter(sun_spectrum, regions=segments)
cutted_sun_spectrum = sun_spectrum[wfilter]



#--- Radial Velocity determination with linelist mask --------------------------
logging.info("Radial velocity determination with linelist mask...")
# - Read atomic data
vald_linelist_file = sve_dir + "/input/linelists/VALD/300_1100nm.lst"
linelist = sve.read_VALD_linelist(vald_linelist_file, minimum_depth=0.0)

###### OPTIONAL:
## - Filter lines that may be affected by telluric lines
#telluric_lines_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
#linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)
#dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
#linelist_telluric = linelist_telluric[dfilter]

#tfilter = sve.create_filter_for_regions_affected_by_tellurics(linelist['wave_peak'], \
                            #linelist_telluric, min_velocity=-30.0, max_velocity=30.0)
#linelist[tfilter]['depth'] = 0.0

xcoord, fluxes, errors = sve.build_velocity_profile(mu_cas_a_spectrum, \
                                            linelist=linelist, lower_velocity_limit=-200.0, \
                                            upper_velocity_limit=200.0, velocity_step=1.0)

models = sve.modelize_velocity_profile(xcoord, fluxes, errors)
best = sve.select_good_velocity_profile_models(models, xcoord, fluxes)
models = models[best]

# Number of models represent the number of components
spectroscopic_nary = str(len(models) >= 2)
# First component:
rv = np.round(models[0].mu(), 2) # km/s


#--- Radial Velocity determination with template -------------------------------
logging.info("Radial velocity determination with template...")
# - Read synthetic template
template = sve.read_spectrum(sve_dir + \
        "/input/spectra/synthetic/Synth_Meszaros_VALD_5777.0_4.44_0.0_1.0.txt.gz")
# - Read observed template
#template = sve.read_spectrum(sve_dir + "/input/spectra/examples/narval_sun.s.gz")

###### OPTIONAL:
## - Read telluric lines and use only the 25% of the deepest ones
#telluric_lines_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
#linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)
#dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
#linelist_telluric = linelist_telluric[dfilter]

## - Filter regions that may be affected by telluric lines
#tfilter = sve.create_filter_for_regions_affected_by_tellurics(template['waveobs'], \
                            #linelist_telluric, min_velocity=-30.0, max_velocity=30.0)
#template['flux'][tfilter] = 0.0


xcoord, fluxes, errors = sve.build_velocity_profile(mu_cas_a_spectrum, \
                                            template=template, lower_velocity_limit=-200.0, \
                                            upper_velocity_limit=200.0, velocity_step=1.0)

models = sve.modelize_velocity_profile(xcoord, fluxes, errors)
best = sve.select_good_velocity_profile_models(models, xcoord, fluxes)
models = models[best]

# Number of models represent the number of components
spectroscopic_nary = str(len(models) >= 2)
# First component:
rv = np.round(models[0].mu(), 2) # km/s


#--- Radial Velocity correction ------------------------------------------------
logging.info("Radial velocity correction...")
mu_cas_a_spectrum = sve.correct_velocity(mu_cas_a_spectrum, rv)


#--- Barycentric Velocity determination ----------------------------------------
logging.info("Barycentric velocity determination...")
# - Telluric
telluric_linelist_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
linelist_telluric = sve.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

xcoord, fluxes, errors = sve.build_velocity_profile(sun_spectrum, \
                                        linelist=linelist_telluric, lower_velocity_limit=-100.0, \
                                        upper_velocity_limit=100.0, velocity_step=0.5)

models = sve.modelize_velocity_profile(xcoord, fluxes, errors, only_one_peak=True)
bv = np.round(models[0].mu(), 2) # km/s


#--- Resolution degradation ----------------------------------------------------
logging.info("Resolution degradation...")
from_resolution = 80000
to_resolution = 40000
convolved_sun_spectrum = sve.convolve_spectrum(sun_spectrum, to_resolution, \
                                                from_resolution=from_resolution)
convolved_mu_cas_a_spectrum = sve.convolve_spectrum(mu_cas_a_spectrum, to_resolution, \
                                                from_resolution=from_resolution)


#--- Resampling and combining --------------------------------------------------
logging.info("Resampling and comibining...")
wavelengths = np.arange(480.0, 670.0, 0.001)
resampled_sun_spectrum = sve.resample_spectrum(sun_spectrum, wavelengths)
resampled_mu_cas_a_spectrum = sve.resample_spectrum(mu_cas_a_spectrum, wavelengths)

# Coadd previously resampled spectra
total_wavelengths = len(resampled_sun_spectrum)
coadded_spectrum = sve.create_spectrum_structure(resampled_sun_spectrum['waveobs'])
coadded_spectrum['flux'] = resampled_sun_spectrum['flux'] + resampled_mu_cas_a_spectrum['flux']
coadded_spectrum['err'] = np.sqrt(np.power(resampled_sun_spectrum['err'],2) + \
                                np.power(resampled_mu_cas_a_spectrum['err'],2))


#--- Fit continuum -------------------------------------------------------------
logging.info("Fit continuum...")

nknots = None # Number of knots will be automatically estimated
model = "Polynomy" # "Splines"

# - Use a fixed value (useful when the spectrum is already normalized)
sun_continuum_model = sve.fit_continuum(sun_spectrum, fixed_value=1.0, model="Fixed value")

# - Consider only continuum regions for finding continuum points to fit
continuum_regions = sve.read_continuum_regions(sve_dir + "/input/regions/fe_lines_continuum.txt")
sun_continuum_model = sve.fit_continuum(sun_spectrum, \
                        continuum_regions=continuum_regions, nknots=nknots,\
                        median_wave_range=0.1, max_wave_range=1.0,
                        model=model)

# - Fit continuum in each segment independently
segments = sve.read_segment_regions(sve_dir + "/input/regions/fe_lines_segments.txt")
sun_continuum_model = sve.fit_continuum(sun_spectrum, \
                        independent_regions=segments, nknots=1,\
                        median_wave_range=0.1, max_wave_range=1.0,
                        model=model)

# - Use the whole spectrum and a polynomial/spline model
sun_continuum_model = sve.fit_continuum(sun_spectrum, nknots=nknots, \
                            median_wave_range=0.1, max_wave_range=1.0, \
                            model=model)


#--- Continuum normalization ---------------------------------------------------
logging.info("Continuum normalization...")
normalized_sun_spectrum = sve.create_spectrum_structure(sun_spectrum['waveobs'])
normalized_sun_spectrum['flux'] = sun_spectrum['flux'] \
                                    / sun_continuum_model(sun_spectrum['waveobs'])
normalized_sun_spectrum['err'] = sun_spectrum['err'] \
                                    / sun_continuum_model(sun_spectrum['waveobs'])


#--- Find continuum regions ----------------------------------------------------
logging.info("Finding continuum regions...")
resolution = 80000
sigma = 0.001
max_continuum_diff = 1.0
fixed_wave_step = 0.05
sun_continuum_regions = sve.find_continuum(sun_spectrum, resolution, \
                                    max_std_continuum = sigma, \
                                    continuum_model = sun_continuum_model, \
                                    max_continuum_diff=max_continuum_diff, \
                                    fixed_wave_step=fixed_wave_step)
sve.write_continuum_regions(sun_continuum_regions, "example_sun_fe_lines_continuum.txt")

# Or limit the search to given segments
segments = sve.read_segment_regions(sve_dir + "/input/regions/fe_lines_segments.txt")
limited_sun_continuum_regions = sve.find_continuum(sun_spectrum, resolution, \
                                        segments=segments, max_std_continuum = sigma, \
                                        continuum_model = sun_continuum_model, \
                                        max_continuum_diff=max_continuum_diff, \
                                        fixed_wave_step=fixed_wave_step)
sve.write_continuum_regions(limited_sun_continuum_regions, \
        "example_limited_sun_continuum_region.txt")

#--- Find linemasks ------------------------------------------------------------
logging.info("Finding line masks...")
sun_continuum_model = sve.fit_continuum(sun_spectrum)
vald_linelist_file = sve_dir + "/input/linelists/VALD/300_1100nm.lst"
chemical_elements_file = sve_dir + "/input/abundances/chemical_elements_symbols.dat"
molecules_file = sve_dir + "/input/abundances/molecular_symbols.dat"
telluric_linelist_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
resolution = 80000
smoothed_sun_spectrum = sve.convolve_spectrum(sun_spectrum, resolution)
min_depth = 0.05
max_depth = 1.00
vel_atomic = 0.00  # km/s
vel_telluric = 17.79  # km/s
sun_linemasks = sve.find_linemasks(sun_spectrum, sun_continuum_model, vald_linelist_file, \
                        chemical_elements_file, molecules_file, telluric_linelist_file, \
                        minimum_depth=min_depth, maximum_depth=max_depth, \
                        smoothed_spectrum=smoothed_sun_spectrum, \
                        discard_gaussian=False, discard_voigt=True, \
                        vel_atomic=vel_atomic, vel_telluric=vel_telluric)
# Exclude lines that have not been successfully cross matched with the atomic data
# because we cannot calculate the chemical abundance (it will crash the corresponding routines)
rejected_by_atomic_line_not_found = (sun_linemasks['VALD_wave_peak'] == 0)
sun_linemasks = sun_linemasks[~rejected_by_atomic_line_not_found]

sve.write_line_regions(sun_linemasks, "example_sun_linemasks.txt")

#--- Barycentric velocity correction -------------------------------------------
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
barycentric_vel = sve.calculate_barycentric_velocity_correction((year, month, day, \
                                hours, minutes, seconds), (ra_hours, ra_minutes, \
                                ra_seconds, dec_degrees, dec_minutes, dec_seconds))
corrected_spectrum = sve.correct_velocity(mu_cas_a_spectrum, barycentric_vel)


#--- Estimate SNR from flux ----------------------------------------------------
logging.info("Estimating SNR...")
num_points = 10
estimated_snr = sve.estimate_snr(sun_spectrum['flux'], num_points=num_points)


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

telluric_lines_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)

# - Filter regions that may be affected by telluric lines
rv = 0.0
min_vel = -30.0
max_vel = +30.0
# Only the 25% of the deepest ones:
dfilter = linelist_telluric['depth'] > np.percentile(linelist_telluric['depth'], 75)
tfilter = sve.create_filter_for_regions_affected_by_tellurics(sun_spectrum['waveobs'], \
                            linelist_telluric[dfilter], min_velocity=-rv+min_vel, \
                            max_velocity=-rv+max_vel)
clean_sun_spectrum = sun_spectrum[tfilter]


#--- Adjust line masks ---------------------------------------------------------
resolution = 80000
smoothed_sun_spectrum = sve.convolve_spectrum(sun_spectrum, resolution)
line_regions = sve.read_line_regions(sve_dir + "/input/regions/fe_lines.txt")
linemasks = sve.adjust_linemasks(smoothed_sun_spectrum, line_regions, margin=0.5)
segments = sve.create_segments_around_lines(linemasks, margin=0.25)


#---Create segments around linemasks -------------------------------------------
line_regions = sve.read_line_regions(sve_dir + "/input/regions/fe_lines.txt")
segments = sve.create_segments_around_lines(line_regions, margin=0.25)


#--- Fit lines -----------------------------------------------------------------
logging.info("Fitting lines...")
vald_linelist_file = sve_dir + "/input/linelists/VALD/300_1100nm.lst"
chemical_elements_file = sve_dir + "/input/abundances/chemical_elements_symbols.dat"
molecules_file = sve_dir + "/input/abundances/molecular_symbols.dat"
telluric_linelist_file = sve_dir + "/input/linelists/telluric/standard_atm_air_model.lst"
vel_atomic = 0.00 # km/s
vel_telluric = 17.79 # km/s
line_regions = sve.read_line_regions(sve_dir + "/input/regions/fe_lines.txt")
linemasks = sve.fit_lines(line_regions, sun_spectrum, sun_continuum_model, vel_atomic, \
                            vel_telluric, vald_linelist_file, chemical_elements_file, \
                            molecules_file, telluric_linelist_file, discard_gaussian=False, \
                            discard_voigt=True)
# Exclude lines that have not been successfully cross matched with the atomic data
# because we cannot calculate the chemical abundance (it will crash the corresponding routines)
rejected_by_atomic_line_not_found = (linemasks['VALD_wave_peak'] == 0)
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
model = "MARCS" # MARCS", "MARCS.GES", "MARCS.APOGEE", "ATLAS9.APOGEE",
                # "ATLAS9.Castelli", "ATLAS9.Kurucz", "ATLAS9.Kirby"
linelist_name = "VALD.300_1100nm"   # "VALD.300_1100nm",
                                    # "VALD_with_molecules.300_1100nm",
                                    # "GES.475_685nm",  "Kurucz.300_1100nm",
                                    # "NIST.300_1100nm", "SPECTRUM.300_1000nm"]
solar_abundances = "Grevesse.2007"  # "Asplund.2005", "Asplund.2009",
                                    # "Grevesse.1998", "Anders.1989"

# Load model atmospheres
modeled_layers_pack = sve.load_modeled_layers_pack(sve_dir + '/input/atmospheres/' \
        + model + '/modeled_layers_pack.dump')
# Load SPECTRUM abundances
abundances_file = sve_dir + "/input/abundances/" + solar_abundances + "/stdatom.dat"
abundances = sve.read_SPECTRUM_abundances(abundances_file)

# Validate parameters
if not sve.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
    msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
            fall out of theatmospheric models."
    print msg

# Prepare atmosphere model
atmosphere_layers = sve.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)

spec_abund, normal_abund, x_over_h, x_over_fe = sve.determine_abundances(atmosphere_layers, \
        teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, verbose=1)

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
#regions = sve.read_segment_regions(sve_dir + "/input/regions/fe_lines_segments.txt")
regions = None
wave_base = 515.0 # Magnesium triplet region
wave_top = 525.0

# Selected model amtosphere, linelist and solar abundances
model = "MARCS" # MARCS", "MARCS.GES", "MARCS.APOGEE", "ATLAS9.APOGEE",
                # "ATLAS9.Castelli", "ATLAS9.Kurucz", "ATLAS9.Kirby"
linelist_name = "VALD.300_1100nm"   # "VALD.300_1100nm",
                                    # "VALD_with_molecules.300_1100nm",
                                    # "GES.475_685nm",  "Kurucz.300_1100nm",
                                    # "NIST.300_1100nm", "SPECTRUM.300_1000nm"]
solar_abundances = "Grevesse.2007"  # "Asplund.2005", "Asplund.2009",
                                    # "Grevesse.1998", "Anders.1989"

# Load model atmospheres
modeled_layers_pack = sve.load_modeled_layers_pack(sve_dir + 'input/atmospheres/' + \
        model + '/modeled_layers_pack.dump')
# Load SPECTRUM linelist
linelist_file = sve_dir + "/input/linelists/SPECTRUM/" + linelist_name.split(".")[0] +\
                    "/" + linelist_name.split(".")[1] + ".lst"
linelist = sve.read_SPECTRUM_linelist(linelist_file)
# Load SPECTRUM abundances
fixed_abundances = None # No fixed abundances
abundances_file = sve_dir + "/input/abundances/" + solar_abundances + "/stdatom.dat"
abundances = sve.read_SPECTRUM_abundances(abundances_file)

# Validate parameters
if not sve.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
    msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
            fall out of theatmospheric models."
    print msg


# Prepare atmosphere model
atmosphere_layers = sve.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH)

# Synthesis
synth_spectrum = sve.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
synth_spectrum['flux'] = sve.generate_spectrum(synth_spectrum['waveobs'], \
        atmosphere_layers, teff, logg, MH, linelist=linelist, abundances=abundances, \
        fixed_abundances=fixed_abundances, microturbulence_vel = microturbulence_vel, \
        macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
        R=resolution, regions=regions, verbose=1)


#--- Adding gaussian noise -----------------------------------------------------
snr = 100
distribution = "poisson" # "gaussian"
synth_spectrum = sve.add_noise(synth_spectrum, snr, distribution)


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
model = "MARCS" # MARCS", "MARCS.GES", "MARCS.APOGEE", "ATLAS9.APOGEE",
                # "ATLAS9.Castelli", "ATLAS9.Kurucz", "ATLAS9.Kirby"
linelist_name = "VALD.300_1100nm"   # "VALD.300_1100nm",
                                    # "VALD_with_molecules.300_1100nm",
                                    # "GES.475_685nm",  "Kurucz.300_1100nm",
                                    # "NIST.300_1100nm", "SPECTRUM.300_1000nm"]
solar_abundances = "Grevesse.2007"  # "Asplund.2005", "Asplund.2009",
                                    # "Grevesse.1998", "Anders.1989"

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
segments = sve.read_segment_regions(sve_dir + "/input/regions/fe_lines_segments.txt")
line_regions = sve.read_line_regions(sve_dir + "/input/regions/fe_lines.txt")

# Load model atmospheres
modeled_layers_pack = sve.load_modeled_layers_pack(sve_dir + 'input/atmospheres/' + \
        model + '/modeled_layers_pack.dump')
# Load SPECTRUM linelist
linelist_file = sve_dir + "/input/linelists/SPECTRUM/" + linelist_name.split(".")[0] +\
                    "/" + linelist_name.split(".")[1] + ".lst"
linelist = sve.read_SPECTRUM_linelist(linelist_file)
# Load SPECTRUM abundances
fixed_abundances = None # No fixed abundances
abundances_file = sve_dir + "/input/abundances/" + solar_abundances + "/stdatom.dat"
abundances = sve.read_SPECTRUM_abundances(abundances_file)

obs_spec, modeled_synth_spectrum, params, errors, free_abundances, status, stats_linemasks = \
        sve.modelize_spectrum(sun_spectrum, sun_continuum_model, \
        modeled_layers_pack, linelist, abundances, free_abundances, initial_teff, \
        initial_logg, initial_MH, initial_vmic, initial_vmac, initial_vsini, \
        initial_limb_darkening_coeff, initial_R, free_params, segments=segments, \
        linemasks=line_regions, max_iterations=max_iterations)


##--- Save spectrum ------------------------------------------------------------
logging.info("Saving spectrum...")
sve.write_spectrum(sun_spectrum, "example_sun.s")
sve.write_spectrum(mu_cas_a_spectrum, "example_mu_cas_a.s")
sve.write_spectrum(synth_spectrum, "example_synth.s")
sve.write_spectrum(modeled_synth_spectrum, "example_modeled_synth.s")


