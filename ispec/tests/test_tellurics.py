import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestTellurics(unittest.TestCase):

    def test_determine_tellurics_shift_with_mask(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Telluric velocity shift determination from spectrum --------------------------
        # - Telluric
        telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

        models, ccf = ispec.cross_correlate_with_mask(star_spectrum, telluric_linelist, \
                                lower_velocity_limit=-100, upper_velocity_limit=100, \
                                velocity_step=0.5, mask_depth=0.01, \
                                fourier = False,
                                only_one_peak = True)

        bv = models[0].mu() # km/s
        bv_err = models[0].emu() # km/s
        self.assertAlmostEqual(bv, 8.139999999999954)
        self.assertAlmostEqual(bv_err, 0.031123202937016047)

    def test_determine_tellurics_shift_with_template(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Telluric velocity shift determination from spectrum --------------------------
        # - Read synthetic template
        template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Tellurics.350_1100nm/template.txt.gz")

        models, ccf = ispec.cross_correlate_with_template(star_spectrum, template, \
                                lower_velocity_limit=-100, upper_velocity_limit=100, \
                                velocity_step=0.5, fourier=False, \
                                only_one_peak = True)

        bv = models[0].mu() # km/s
        bv_err = models[0].emu() # km/s
        self.assertAlmostEqual(bv, 8.169999999999954)
        self.assertAlmostEqual(bv_err, 0.6605884013888695)

    def test_clean_telluric_regions(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #--- Telluric velocity shift determination from spectrum --------------------------
        # - Telluric
        telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

        models, ccf = ispec.cross_correlate_with_mask(star_spectrum, telluric_linelist, \
                                lower_velocity_limit=-100, upper_velocity_limit=100, \
                                velocity_step=0.5, mask_depth=0.01, \
                                fourier = False,
                                only_one_peak = True)

        bv = np.round(models[0].mu(), 2) # km/s
        bv_err = np.round(models[0].emu(), 2) # km/s

        #--- Clean regions that may be affected by tellurics ---------------------------
        telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
        telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

        # - Filter regions that may be affected by telluric lines
        #bv = 0.0
        min_vel = -30.0
        max_vel = +30.0
        # Only the 25% of the deepest ones:
        dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
        tfilter = ispec.create_filter_for_regions_affected_by_tellurics(star_spectrum['waveobs'], \
                                    telluric_linelist[dfilter], min_velocity=-bv+min_vel, \
                                    max_velocity=-bv+max_vel)
        clean_star_spectrum = star_spectrum[~tfilter]
        self.assertEqual(len(np.where(tfilter)[0]), 686)
        self.assertEqual(np.where(tfilter)[0][0], 45002)
        self.assertEqual(np.where(tfilter)[0][-1], 46051)
