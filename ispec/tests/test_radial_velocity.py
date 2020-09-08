import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestRadialVelocity(unittest.TestCase):

    def test_determine_radial_velocity_with_mask(self):
        mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
        #--- Radial Velocity determination with linelist mask --------------------------
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
        ccf_mask = ispec.read_cross_correlation_mask(mask_file)

        models, ccf = ispec.cross_correlate_with_mask(mu_cas_spectrum, ccf_mask, \
                                lower_velocity_limit=-200, upper_velocity_limit=200, \
                                velocity_step=1.0, mask_depth=0.01, \
                                fourier=False)

        # Number of models represent the number of components
        components = len(models)
        # First component:
        rv = models[0].mu() # km/s
        rv_err = models[0].emu() # km/s
        self.assertEqual(components, 1)
        self.assertAlmostEqual(rv, -96.43)
        self.assertAlmostEqual(rv_err, 0.03736164048308909)

    def determine_radial_velocity_with_template(self):
        mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
        #--- Radial Velocity determination with template -------------------------------
        # - Read synthetic template
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Arcturus.372_926nm/template.txt.gz")
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Sun.372_926nm/template.txt.gz")
        template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
        #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

        models, ccf = ispec.cross_correlate_with_template(mu_cas_spectrum, template, \
                                lower_velocity_limit=-200, upper_velocity_limit=200, \
                                velocity_step=1.0, fourier=False)

        # Number of models represent the number of components
        components = len(models)
        # First component:
        rv = models[0].mu() # km/s
        rv_err = models[0].emu() # km/s
        self.assertEqual(components, 1)
        self.assertAlmostEqual(rv, -96.43)
        self.assertAlmostEqual(rv_err, 0.03736164048308909)

    def test_correct_radial_velocity(self):
        mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
        #--- Radial Velocity correction ------------------------------------------------
        rv = -96.40 # km/s
        mu_cas_spectrum = ispec.correct_velocity(mu_cas_spectrum, rv)
        self.assertAlmostEqual(mu_cas_spectrum['waveobs'][0], 480.15547196)
        self.assertAlmostEqual(mu_cas_spectrum['flux'][0], 0.19076)
        self.assertAlmostEqual(mu_cas_spectrum['err'][0], 0.00095993)

