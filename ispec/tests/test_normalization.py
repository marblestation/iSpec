import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestNormalization(unittest.TestCase):

    def test_normalize_spectrum_using_continuum_regions(self):
        """
        Consider only continuum regions for the fit, strategy 'median+max'
        """
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

        continuum_regions = ispec.read_continuum_regions(ispec_dir + "/input/regions/fe_lines_continuum.txt")
        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                continuum_regions=continuum_regions, nknots=nknots, degree=degree, \
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

        #--- Continuum normalization ---------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")

        np.testing.assert_equal(star_spectrum['waveobs'], normalized_star_spectrum['waveobs'])
        self.assertAlmostEqual(star_spectrum['flux'][0], 0.67077)
        self.assertAlmostEqual(star_spectrum['flux'][-1], 2.2169)
        self.assertAlmostEqual(star_spectrum['err'][0], 0.0021259)
        self.assertAlmostEqual(star_spectrum['err'][-1], 0.0043878)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][0], 0.7990890725070141)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][-1], 0.9887281323707806)
        self.assertAlmostEqual(normalized_star_spectrum['err'][0], 0.0025325871151701197)
        self.assertAlmostEqual(normalized_star_spectrum['err'][-1], 0.001956940457042046)
        self.assertTrue(np.all(star_continuum_model(star_spectrum['waveobs']) == 1))

    def test_normalize_spectrum_in_segments(self):
        """
        Fit continuum in each segment independently, strategy 'median+max'
        """
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")

        #--- Continuum fit -------------------------------------------------------------
        model = "Splines" # "Polynomy"
        degree = 2
        nknots = 1
        from_resolution = 80000

        # Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
        order='median+max'
        median_wave_range=0.05
        max_wave_range=1.0

        segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                independent_regions=segments, nknots=nknots, degree=degree,\
                                median_wave_range=median_wave_range, \
                                max_wave_range=max_wave_range, \
                                model=model, order=order, \
                                automatic_strong_line_detection=True, \
                                strong_line_probability=0.5, \
                                use_errors_for_fitting=True)

        #--- Continuum normalization ---------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        np.testing.assert_equal(star_spectrum['waveobs'], normalized_star_spectrum['waveobs'])
        self.assertAlmostEqual(star_spectrum['flux'][0], 0.67077)
        self.assertAlmostEqual(star_spectrum['flux'][-1], 2.2169)
        self.assertAlmostEqual(star_spectrum['err'][0], 0.0021259)
        self.assertAlmostEqual(star_spectrum['err'][-1], 0.0043878)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][0], 1.0)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][-1], 1.0)
        self.assertAlmostEqual(normalized_star_spectrum['err'][0], 0.0)
        self.assertAlmostEqual(normalized_star_spectrum['err'][-1], 0.0)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][11], 0.9928764307623508)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][58288], 0.960387654522309)
        self.assertAlmostEqual(normalized_star_spectrum['err'][11], 0.0028468129511384785)
        self.assertAlmostEqual(normalized_star_spectrum['err'][58288], 0.0019771951824093786)
        self.assertEqual(len(np.where(normalized_star_spectrum['flux'] != 1.)[0]), 33722)
        self.assertTrue(np.all(star_continuum_model(star_spectrum['waveobs']) == 1))

    def test_normalize_whole_spectrum_with_template(self):
        """
        Use a template to normalize the whole spectrum
        """
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        synth_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

        #--- Continuum fit -------------------------------------------------------------
        model = "Template"
        nknots = None # Automatic: 1 spline every 5 nm (in this case, used to apply a gaussian filter)
        from_resolution = 80000
        median_wave_range=5.0

        #strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
        strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/relevant/relevant_line_masks.txt")
        #strong_lines = None
        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                    ignore=strong_lines, \
                                    nknots=nknots, \
                                    median_wave_range=median_wave_range, \
                                    model=model, \
                                    template=synth_spectrum)

        #--- Continuum normalization ---------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        np.testing.assert_equal(star_spectrum['waveobs'], normalized_star_spectrum['waveobs'])
        self.assertAlmostEqual(star_spectrum['flux'][0], 0.67077)
        self.assertAlmostEqual(star_spectrum['flux'][-1], 2.2169)
        self.assertAlmostEqual(star_spectrum['err'][0], 0.0021259)
        self.assertAlmostEqual(star_spectrum['err'][-1], 0.0043878)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][0], 0.7964201448878107)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][-1], 0.997168332612696)
        self.assertAlmostEqual(normalized_star_spectrum['err'][0], 0.002524128368914824)
        self.assertAlmostEqual(normalized_star_spectrum['err'][-1], 0.0019736457259407225)
        self.assertTrue(np.all(star_continuum_model(star_spectrum['waveobs']) == 1))

    def test_normalize_whole_spectrum(self):
        """
        Use the whole spectrum, strategy 'median+max'
        """
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

        #--- Continuum normalization ---------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        np.testing.assert_equal(star_spectrum['waveobs'], normalized_star_spectrum['waveobs'])
        self.assertAlmostEqual(star_spectrum['flux'][0], 0.67077)
        self.assertAlmostEqual(star_spectrum['flux'][-1], 2.2169)
        self.assertAlmostEqual(star_spectrum['err'][0], 0.0021259)
        self.assertAlmostEqual(star_spectrum['err'][-1], 0.0043878)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][0], 0.8052457389161867)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][-1], 0.9895767577073404)
        self.assertAlmostEqual(normalized_star_spectrum['err'][0], 0.0025520997008839415)
        self.assertAlmostEqual(normalized_star_spectrum['err'][-1], 0.0019586200989978207)
        self.assertTrue(np.all(star_continuum_model(star_spectrum['waveobs']) == 1))

    def test_normalize_whole_spectrum_ignoring_prefixed_strong_lines(self):
        """
        Use the whole spectrum but ignoring some strong lines, strategy 'median+max'
        """
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

        strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
        star_continuum_model = ispec.fit_continuum(star_spectrum, from_resolution=from_resolution, \
                                    ignore=strong_lines, \
                                    nknots=nknots, degree=degree, \
                                    median_wave_range=median_wave_range, \
                                    max_wave_range=max_wave_range, \
                                    model=model, order=order, \
                                    automatic_strong_line_detection=True, \
                                    strong_line_probability=0.5, \
                                    use_errors_for_fitting=True)

        #--- Continuum normalization ---------------------------------------------------
        normalized_star_spectrum = ispec.normalize_spectrum(star_spectrum, star_continuum_model, consider_continuum_errors=False)
        # Use a fixed value because the spectrum is already normalized
        star_continuum_model = ispec.fit_continuum(star_spectrum, fixed_value=1.0, model="Fixed value")
        np.testing.assert_equal(star_spectrum['waveobs'], normalized_star_spectrum['waveobs'])
        self.assertAlmostEqual(star_spectrum['flux'][0], 0.67077)
        self.assertAlmostEqual(star_spectrum['flux'][-1], 2.2169)
        self.assertAlmostEqual(star_spectrum['err'][0], 0.0021259)
        self.assertAlmostEqual(star_spectrum['err'][-1], 0.0043878)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][0], 0.8005567874185445)
        self.assertAlmostEqual(normalized_star_spectrum['flux'][-1], 0.9874528805693361)
        self.assertAlmostEqual(normalized_star_spectrum['err'][0], 0.0025372388067043607)
        self.assertAlmostEqual(normalized_star_spectrum['err'][-1], 0.0019544164145257493)
        self.assertTrue(np.all(star_continuum_model(star_spectrum['waveobs']) == 1))
