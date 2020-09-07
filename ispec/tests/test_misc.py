import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestMisc(unittest.TestCase):

    def test_calculate_barycentric_velocity(self):
        mu_cas_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
        #--- Barycentric velocity correction from observation date/coordinates ---------
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
        self.assertAlmostEquals(barycentric_vel, -10.32)
