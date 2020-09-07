import os
import sys
import unittest
import numpy as np
import tempfile

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestRegions(unittest.TestCase):

    def test_find_continuum_regions(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
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
        #--- Find continuum regions ----------------------------------------------------
        resolution = 80000
        sigma = 0.001
        max_continuum_diff = 1.0
        fixed_wave_step = 0.05
        star_continuum_regions = ispec.find_continuum(star_spectrum, resolution, \
                                            max_std_continuum = sigma, \
                                            continuum_model = star_continuum_model, \
                                            max_continuum_diff=max_continuum_diff, \
                                            fixed_wave_step=fixed_wave_step)
        #ispec.write_continuum_regions(star_continuum_regions, "example_star_fe_lines_continuum.txt")
        self.assertEquals(len(star_continuum_regions), 2)
        self.assertAlmostEquals(star_continuum_regions['wave_base'][0], 528.6665701)
        self.assertAlmostEquals(star_continuum_regions['wave_top'][0], 528.6915701)
        self.assertAlmostEquals(star_continuum_regions['wave_base'][1], 608.4415701)
        self.assertAlmostEquals(star_continuum_regions['wave_top'][1], 608.4665701)


    def test_find_continuum_regions_in_segments(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
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
        #--- Find continuum regions in segments ----------------------------------------
        resolution = 80000
        sigma = 0.005
        max_continuum_diff = 1.0
        fixed_wave_step = 0.05
        # Limit the search to given segments
        segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        limited_star_continuum_regions = ispec.find_continuum(star_spectrum, resolution, \
                                                segments=segments, max_std_continuum = sigma, \
                                                continuum_model = star_continuum_model, \
                                                max_continuum_diff=max_continuum_diff, \
                                                fixed_wave_step=fixed_wave_step)
        #ispec.write_continuum_regions(limited_star_continuum_regions, \
                #"example_limited_star_continuum_region.txt")
        self.assertEquals(len(limited_star_continuum_regions), 504)


    def test_find_linemasks(self):
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
        #--- Find linemasks ------------------------------------------------------------
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

        telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"

        # Read
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun
        #telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)
        #vel_telluric = 17.79  # km/s
        #telluric_linelist = None
        #vel_telluric = None

        resolution = 80000
        smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, resolution)
        min_depth = 0.05
        max_depth = 1.00
        star_linemasks = ispec.find_linemasks(star_spectrum, star_continuum_model, \
                                atomic_linelist=atomic_linelist, \
                                max_atomic_wave_diff = 0.005, \
                                telluric_linelist=telluric_linelist, \
                                vel_telluric=vel_telluric, \
                                minimum_depth=min_depth, maximum_depth=max_depth, \
                                smoothed_spectrum=smoothed_star_spectrum, \
                                check_derivatives=False, \
                                discard_gaussian=False, discard_voigt=True, \
                                closest_match=False)
        # Exclude lines that have not been successfully cross matched with the atomic data
        # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
        rejected_by_atomic_line_not_found = (star_linemasks['wave_nm'] == 0)
        star_linemasks = star_linemasks[~rejected_by_atomic_line_not_found]

        # Exclude lines with EW equal to zero
        rejected_by_zero_ew = (star_linemasks['ew'] == 0)
        star_linemasks = star_linemasks[~rejected_by_zero_ew]

        # Select only iron lines
        iron = star_linemasks['element'] == "Fe 1"
        iron = np.logical_or(iron, star_linemasks['element'] == "Fe 2")
        iron_star_linemasks = star_linemasks[iron]

        self.assertEquals(len(star_linemasks), 1731)
        self.assertEquals(len(iron_star_linemasks), 883)
        self.assertAlmostEquals(iron_star_linemasks['ew'][0], 15.82970157775808)
        self.assertAlmostEquals(iron_star_linemasks['ew'][-1], 15.545342323635188)
        self.assertAlmostEquals(iron_star_linemasks['ew_err'][0], 1.1187509030236669)
        self.assertAlmostEquals(iron_star_linemasks['ew_err'][-1], 1.3221015426240295)

        tmp_filename = tempfile.mktemp()
        # Write regions with masks limits, cross-matched atomic data and fit data
        ispec.write_line_regions(star_linemasks, tmp_filename, extended=True)
        recover_star_linemasks = ispec.read_line_regions(tmp_filename)
        os.remove(tmp_filename)
        # Write regions with masks limits and cross-matched atomic data (fit data fields are zeroed)
        zeroed_star_linemasks = ispec.reset_fitted_data_fields(star_linemasks)

        np.testing.assert_almost_equal(star_linemasks['wave_nm'], recover_star_linemasks['wave_nm'])
        np.testing.assert_almost_equal(star_linemasks['ew'], recover_star_linemasks['ew'], decimal=1)
        np.testing.assert_almost_equal(star_linemasks['wave_base'], recover_star_linemasks['wave_base'], decimal=4)
        np.testing.assert_almost_equal(star_linemasks['wave_top'], recover_star_linemasks['wave_top'], decimal=4)
        np.testing.assert_almost_equal(star_linemasks['wave_peak'], recover_star_linemasks['wave_peak'], decimal=4)
        self.assertTrue(np.all(zeroed_star_linemasks['ew'] == 0))

    def test_adjust_line_masks(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Adjust line masks ---------------------------------------------------------
        resolution = 80000
        smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, resolution)
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
        linemasks = ispec.adjust_linemasks(smoothed_star_spectrum, line_regions, max_margin=0.5)
        self.assertAlmostEquals(line_regions['wave_base'][0], 480.26937)
        self.assertAlmostEquals(line_regions['wave_top'][0], 480.28771)
        self.assertEqual(len(linemasks), 315)
        np.testing.assert_equal(line_regions['wave_peak'], linemasks['wave_peak'])
        self.assertAlmostEquals(linemasks['wave_base'][0], 480.269365109)
        self.assertAlmostEquals(linemasks['wave_top'][0], 480.319964199)

    def test_create_segments_around_linemasks(self):
        #---Create segments around linemasks -------------------------------------------
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/fe_lines.txt")
        segments = ispec.create_segments_around_lines(line_regions, margin=0.25)
        self.assertEquals(len(segments), 132)
        self.assertAlmostEquals(segments['wave_base'][0], 480.01937)
        self.assertAlmostEquals(segments['wave_top'][0], 481.08295)
