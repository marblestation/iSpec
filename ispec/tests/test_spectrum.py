import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestSpectrum(unittest.TestCase):

    def test_convert_air_to_vacuum(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Converting wavelengths from air to vacuum and viceversa -------------------
        star_spectrum_vacuum = ispec.air_to_vacuum(star_spectrum)
        star_spectrum_air = ispec.vacuum_to_air(star_spectrum_vacuum)
        np.testing.assert_equal(star_spectrum['flux'], star_spectrum_vacuum['flux'])
        np.testing.assert_equal(star_spectrum['err'], star_spectrum_vacuum['err'])
        np.testing.assert_equal(star_spectrum['flux'], star_spectrum_air['flux'])
        np.testing.assert_equal(star_spectrum['err'], star_spectrum_air['err'])
        np.testing.assert_almost_equal(star_spectrum['waveobs'], star_spectrum_air['waveobs'])
        self.assertAlmostEqual(star_spectrum_vacuum['waveobs'][0], 480.12574305216896)
        self.assertAlmostEqual(star_spectrum_vacuum['waveobs'][-1], 680.175349351457)

    def test_cut_spectrum_from_segments(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Cut -----------------------------------------------------------------------
        # Keep only points inside a list of segments
        segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        wfilter = ispec.create_wavelength_filter(star_spectrum, regions=segments)
        cutted_star_spectrum = star_spectrum[wfilter]
        self.assertEqual(len(segments), 132)
        self.assertAlmostEqual(segments['wave_base'][0], 480.01937)
        self.assertAlmostEqual(segments['wave_top'][0], 481.08295)
        self.assertEqual(len(cutted_star_spectrum), 33722)
        self.assertAlmostEqual(cutted_star_spectrum['waveobs'][0], 480.02156956)
        self.assertAlmostEqual(cutted_star_spectrum['flux'][0], 0.81984)
        self.assertAlmostEqual(cutted_star_spectrum['err'][0], 0.0023458)

    def test_degrade_resolution(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Resolution degradation ----------------------------------------------------
        from_resolution = 80000
        to_resolution = 47000
        convolved_star_spectrum = ispec.convolve_spectrum(star_spectrum, to_resolution, \
                                                        from_resolution=from_resolution)
        np.testing.assert_equal(star_spectrum['waveobs'], convolved_star_spectrum['waveobs'])
        self.assertAlmostEqual(convolved_star_spectrum['flux'][0], 0.7109576121207318)
        self.assertAlmostEqual(convolved_star_spectrum['flux'][-1], 2.2237382541838007)
        self.assertAlmostEqual(convolved_star_spectrum['err'][0], 0.0021974272902411723)
        self.assertAlmostEqual(convolved_star_spectrum['err'][-1], 0.004453852796214372)

    def test_smooth_spectrum(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Smoothing spectrum (resolution will be affected) --------------------------
        resolution = 80000
        smoothed_star_spectrum = ispec.convolve_spectrum(star_spectrum, resolution)
        np.testing.assert_equal(star_spectrum['waveobs'], smoothed_star_spectrum['waveobs'])
        self.assertAlmostEqual(smoothed_star_spectrum['flux'][0], 0.6973069596669491)
        self.assertAlmostEqual(smoothed_star_spectrum['flux'][-1], 2.2207560753191666)
        self.assertAlmostEqual(smoothed_star_spectrum['err'][0], 0.002174423029777254)
        self.assertAlmostEqual(smoothed_star_spectrum['err'][-1], 0.004451564365118501)

    def test_resample_spectrum(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Resampling  --------------------------------------------------------------
        wavelengths = np.arange(480.0, 680.0, 0.001)
        resampled_star_spectrum = ispec.resample_spectrum(star_spectrum, wavelengths, method="linear", zero_edges=True)
        #resampled_star_spectrum = ispec.resample_spectrum(star_spectrum, wavelengths, method="bessel", zero_edges=True)
        self.assertEqual(len(resampled_star_spectrum), 200000)
        self.assertAlmostEqual(resampled_star_spectrum['waveobs'][0], 480.)
        self.assertAlmostEqual(resampled_star_spectrum['waveobs'][-1], 680.-0.001)
        self.assertAlmostEqual(resampled_star_spectrum['flux'][0], 0.7957469310344828)
        self.assertAlmostEqual(resampled_star_spectrum['flux'][-1], 0.0)
        self.assertAlmostEqual(resampled_star_spectrum['err'][0], 0.002313683448275862)
        self.assertAlmostEqual(resampled_star_spectrum['err'][-1], 0.0)

    def test_filter_cosmic_rays(self):
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
        #--- Filtering cosmic rays -----------------------------------------------------
        # Spectrum should be already normalized
        cosmics = ispec.create_filter_cosmic_rays(star_spectrum, star_continuum_model, \
                                                resampling_wave_step=0.001, window_size=15, \
                                                variation_limit=0.01)
        clean_star_spectrum = star_spectrum[~cosmics]
        np.testing.assert_equal(np.where(cosmics)[0], np.array([  730,  1472,  2626,  2635,  2703,  2704,  5319,  5773,  5782, 6075,  6500,  6610, 11814, 11846, 11883, 12050, 12104, 16562, 16862, 25310, 25311, 25948, 29402, 32355, 34342, 46984, 49469, 57858]))

    def test_add_noise_to_spectrum(self):
        """
        Add noise to an spectrum (ideally to a synthetic one) based on a given SNR.
        """
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        np.random.seed(42)
        #--- Adding poisson noise -----------------------------------------------------
        snr = 100
        distribution = "poisson" # "gaussian"
        noisy_star_spectrum = ispec.add_noise(star_spectrum, snr, distribution)
        np.testing.assert_almost_equal(noisy_star_spectrum['flux'][:10], np.array([0.6678, 0.7459, 0.776 , 0.7796, 0.802 , 0.7953, 0.7383, 0.7164, 0.7381, 0.7876]))
        np.testing.assert_equal(star_spectrum['waveobs'], noisy_star_spectrum['waveobs'])

    def test_generate_new_random_realizations_from_spectrum(self):
        """
        Considering fluxes as mean values and errors as standard deviation, generate
        N new random spectra.
        """
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        np.random.seed(42)

        number = 10
        random_realization_spectra = ispec.random_realizations(star_spectrum, number, distribution="poisson")
        self.assertEqual(len(random_realization_spectra), number)
        np.testing.assert_almost_equal(random_realization_spectra[0]['flux'][:10], np.array([0.67014417, 0.74121144, 0.78384863, 0.79379749, 0.80003082, 0.78786734, 0.74521218, 0.72019732, 0.73934592, 0.78538714]))
        for i in range(number):
            np.testing.assert_equal(random_realization_spectra[i]['waveobs'], star_spectrum['waveobs'])
            np.testing.assert_equal(random_realization_spectra[i]['err'], star_spectrum['err'])

    def test_estimate_snr_from_flux(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        ## WARNING: To compare SNR estimation between different spectra, they should
        ##          be homogeneously sampled (consider a uniform re-sampling)
        #--- Estimate SNR from flux ----------------------------------------------------
        num_points = 10
        estimated_snr = ispec.estimate_snr(star_spectrum['flux'], num_points=num_points)
        self.assertAlmostEqual(estimated_snr, 139.92497450174938)

