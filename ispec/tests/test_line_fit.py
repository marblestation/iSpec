import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestLineFit(unittest.TestCase):

    def test_fit_lines_determine_ew_and_crossmatch_with_atomic_data(self):
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
        from_resolution = 80000

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
        #--- Fit lines -----------------------------------------------------------------
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"


        # Read
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

        #telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
        #vel_telluric = 17.79 # km/s
        #telluric_linelist = None
        #vel_telluric = None

        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/moog_synth_good_for_params_all.txt")
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/width_synth_good_for_params_all.txt")
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/moog_synth_good_for_params_all.txt")
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/width_synth_good_for_params_all.txt")
        line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)

        linemasks = ispec.fit_lines(line_regions, normalized_star_spectrum, star_continuum_model, \
                                    atomic_linelist = atomic_linelist, \
                                    #max_atomic_wave_diff = 0.005, \
                                    max_atomic_wave_diff = 0.00, \
                                    telluric_linelist = telluric_linelist, \
                                    smoothed_spectrum = None, \
                                    check_derivatives = False, \
                                    vel_telluric = vel_telluric, discard_gaussian=False, \
                                    discard_voigt=True, \
                                    free_mu=True, crossmatch_with_mu=False, closest_match=False)
        # Discard lines that are not cross matched with the same original element stored in the note
        linemasks = linemasks[linemasks['element'] == line_regions['note']]

        # Exclude lines that have not been successfully cross matched with the atomic data
        # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
        rejected_by_atomic_line_not_found = (linemasks['wave_nm'] == 0)
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
            #linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="0.995", tmp_dir=None, verbose=0)
            #linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="3;5764,5766,6047,6052,6068,6076", tmp_dir=None, verbose=0)
            snr = 50
            linemasks = ispec.update_ew_with_ares(normalized_star_spectrum, linemasks, rejt="%s" % (snr), tmp_dir=None, verbose=0)

        self.assertEqual(len(linemasks), 281)
        self.assertAlmostEqual(linemasks['ew'][0], 68.48284589709996)
        self.assertAlmostEqual(linemasks['ew'][-3], 46.17583047097995)
        self.assertEqual(linemasks['element'][0], 'Fe 1')
        self.assertEqual(linemasks['element'][-3], 'Si 1')
        self.assertAlmostEqual(linemasks['loggf'][0], -1.028)
        self.assertAlmostEqual(linemasks['loggf'][-3], -1.062)

    def test_fit_lines_already_crossmatched_with_atomic_data_and_determine_ew(self):
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

        #--- Read lines with atomic data ------------------------------------------------
        line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/moog_synth_good_for_params_all_extended.txt")
        #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/width_synth_good_for_params_all_extended.txt")
        #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/moog_synth_good_for_params_all_extended.txt")
        #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/width_synth_good_for_params_all_extended.txt")

        smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, 2*to_resolution)
        line_regions_with_atomic_data = ispec.adjust_linemasks(smoothed_star_spectrum, line_regions_with_atomic_data, max_margin=0.5)

        #telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
        #vel_telluric = 17.79 # km/s
        #telluric_linelist = None
        #vel_telluric = None

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

        self.assertEqual(len(linemasks), 281)
        self.assertAlmostEqual(linemasks['ew'][0], 68.62459244466727)
        self.assertAlmostEqual(linemasks['ew'][-3], 46.23135207295078)
        self.assertEqual(linemasks['element'][0], 'Fe 1')
        self.assertEqual(linemasks['element'][-3], 'Si 1')
        self.assertAlmostEqual(linemasks['loggf'][0], -1.028)
        self.assertAlmostEqual(linemasks['loggf'][-3], -1.062)
