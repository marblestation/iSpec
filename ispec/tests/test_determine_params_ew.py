import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestDetermineParamsEW(unittest.TestCase):

    def test_determine_astrophysical_parameters_from_ew(self):
        code = "moog"
        use_ares = False
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Radial Velocity determination with template -------------------------------
        # - Read synthetic template
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Arcturus.372_926nm/template.txt.gz")
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Sun.372_926nm/template.txt.gz")
        template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

        models, ccf = ispec.cross_correlate_with_template(star_spectrum, template, \
                                lower_velocity_limit=-200, upper_velocity_limit=200, \
                                velocity_step=1.0, fourier=False)

        # Number of models represent the number of components
        components = len(models)
        # First component:
        rv = np.round(models[0].mu(), 2) # km/s
        rv_err = np.round(models[0].emu(), 2) # km/s
        #--- Radial Velocity correction ------------------------------------------------
        star_spectrum = ispec.correct_velocity(star_spectrum, rv)
        #--- Telluric velocity shift determination from spectrum --------------------------
        # - Telluric
        telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

        models, ccf = ispec.cross_correlate_with_mask(star_spectrum, telluric_linelist, \
                                lower_velocity_limit=-100, upper_velocity_limit=100, \
                                velocity_step=0.5, mask_depth=0.01, \
                                fourier = False,
                                only_one_peak = True)

        vel_telluric = np.round(models[0].mu(), 2) # km/s
        vel_telluric_err = np.round(models[0].emu(), 2) # km/s
        #--- Resolution degradation ----------------------------------------------------
        # NOTE: The line selection was built based on a solar spectrum with R ~ 47,000 and GES/VALD atomic linelist.
        from_resolution = 80000
        to_resolution = 47000
        star_spectrum = ispec.convolve_spectrum(star_spectrum, to_resolution, from_resolution)
        #--- Continuum fit -------------------------------------------------------------
        model = "Splines" # "Polynomy"
        degree = 2
        nknots = None # Automatic: 1 spline every 5 nm
        #from_resolution = 80000
        from_resolution = to_resolution

        # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
        order='median+max'
        median_wave_range=0.05
        max_wave_range=1.0

        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                    nknots=nknots, degree=degree, \
                                    median_wave_range=median_wave_range, \
                                    max_wave_range=max_wave_range, \
                                    model=model, order=order, \
                                    automatic_strong_line_detection=True, \
                                    strong_line_probability=0.5, \
                                    use_errors_for_fitting=True)
        #--- Normalize -------------------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")

        #telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
        #vel_telluric = 17.79 # km/s
        #telluric_linelist = None
        #vel_telluric = None

        #--- Read lines and adjust them ------------------------------------------------
        if code in ['width', 'moog']:
            line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_ew_ispec_good_for_params_all_extended.txt".format(code))
            #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_ew_ispec_good_for_params_all_extended.txt".format(code))
        else:
            line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
            #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))

        # Select only iron lines
        line_regions_with_atomic_data = line_regions_with_atomic_data[np.logical_or(line_regions_with_atomic_data['note'] == "Fe 1", line_regions_with_atomic_data['note'] == "Fe 2")]

        smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, 2*to_resolution)
        line_regions_with_atomic_data = ispec.adjust_linemasks(smoothed_star_spectrum, line_regions_with_atomic_data, max_margin=0.5)

        #--- Fit the lines but do NOT cross-match with any atomic linelist since they already have that information
        linemasks = ispec.fit_lines(line_regions_with_atomic_data, normalized_star_spectrum, star_continuum_model, \
                                    atomic_linelist = None, \
                                    max_atomic_wave_diff = 0.005, \
                                    telluric_linelist = telluric_linelist, \
                                    check_derivatives = False, \
                                    vel_telluric = vel_telluric, discard_gaussian=False, \
                                    smoothed_spectrum=None, \
                                    discard_voigt=True, \
                                    free_mu=True, crossmatch_with_mu=False, closest_match=False)

        # Discard bad masks
        flux_peak = normalized_star_spectrum['flux'][linemasks['peak']]
        flux_base = normalized_star_spectrum['flux'][linemasks['base']]
        flux_top = normalized_star_spectrum['flux'][linemasks['top']]
        bad_mask = np.logical_or(linemasks['wave_peak'] <= linemasks['wave_base'], linemasks['wave_peak'] >= linemasks['wave_top'])
        bad_mask = np.logical_or(bad_mask, flux_peak >= flux_base)
        bad_mask = np.logical_or(bad_mask, flux_peak >= flux_top)
        linemasks = linemasks[~bad_mask]

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
            #linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="0.995", tmp_dir=None, verbose=0)
            #linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="3;5764,5766,6047,6052,6068,6076", tmp_dir=None, verbose=0)
            snr = 50
            linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="%s" % (snr), tmp_dir=None, verbose=0)


        #--- Model spectra from EW --------------------------------------------------
        # Parameters
        initial_teff = 5777.0
        initial_logg = 4.44
        initial_MH = 0.00
        initial_alpha = 0.00
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        max_iterations = 10

        # Selected model amtosphere, linelist and solar abundances
        #model = ispec_dir + "/input/atmospheres/MARCS/"
        model = ispec_dir + "/input/atmospheres/MARCS.GES/"
        #model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/"
        #model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/"
        #model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/"
        #model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/"
        #model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/"

        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

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
        if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha}):
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
        efilter = np.logical_and(efilter, linemasks['lower_state_eV'] <= 5.0)
        efilter = np.logical_and(efilter, linemasks['lower_state_eV'] >= 0.5)
        ## Filter also bad fits
        efilter = np.logical_and(efilter, linemasks['rms'] < 1.00)
        # no flux
        noflux = normalized_star_spectrum['flux'][linemasks['peak']] < 1.0e-10
        efilter = np.logical_and(efilter, np.logical_not(noflux))
        unfitted = linemasks['fwhm'] == 0
        efilter = np.logical_and(efilter, np.logical_not(unfitted))

        results = ispec.model_spectrum_from_ew(linemasks[efilter], modeled_layers_pack, \
                            solar_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, \
                            free_params=["teff", "logg", "vmic"], \
                            adjust_model_metalicity=True, \
                            max_iterations=max_iterations, \
                            enhance_abundances=True, \
                            #outliers_detection = "robust", \
                            #outliers_weight_limit = 0.90, \
                            outliers_detection = "sigma_clipping", \
                            #sigma_level = 3, \
                            tmp_dir = None, \
                            code=code)
        params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks = results

        expected_params = {
                            'MH': 0.02800000000000047,
                            'alpha': 0.0,
                            'logg': 4.394728420202205,
                            'teff': 5816.815035823323,
                            'vmic': 1.15562084457677
                           }
        for k, v in expected_params.items():
            self.assertAlmostEquals(params[k], v)
        expected_errors = {
                            'MH': 0.06311918023966583,
                            'alpha': 0.0,
                            'logg': 0.08667507551783957,
                            'teff': 60.97370070760989,
                            'vmic': 0.03611148017625298
                          }
        for k, v in expected_errors.items():
            self.assertAlmostEquals(errors[k], v)

