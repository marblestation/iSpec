import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestReader(unittest.TestCase):

    def test_calculate_theoretical_ew_and_depth(self):
        code = "spectrum"
        #--- Calculate theoretical equivalent widths and depths for a linelist ---------
        # Parameters
        teff = 5777.0
        logg = 4.44
        MH = 0.00
        alpha = 0.00
        microturbulence_vel = 1.0

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
        atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv"
        #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv6_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

        # Load chemical information and linelist
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file)
        atomic_linelist = atomic_linelist[:100] # Select only the first 100 lines (just to reduce execution time, don't do it in a real analysis)

        isotopes = ispec.read_isotope_data(isotope_file)

        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)

        # Load SPECTRUM abundances
        solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

        # Validate parameters
        if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                    fall out of theatmospheric models."
            print(msg)

        # Prepare atmosphere model
        atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha})

        # Synthesis
        #output_wave, output_code, output_ew, output_depth = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
        new_atomic_linelist = ispec.calculate_theoretical_ew_and_depth(atmosphere_layers, \
                teff, logg, MH, alpha, \
                atomic_linelist[:10], isotopes, solar_abundances, microturbulence_vel=microturbulence_vel, \
                verbose=1, gui_queue=None, timeout=900)
        #ispec.write_atomic_linelist(new_atomic_linelist, linelist_filename="example_linelist.txt")

        np.testing.assert_almost_equal(new_atomic_linelist['theoretical_ew'][:10], np.array([2.500e-01, 1.200e+00, 3.000e-02, 2.000e-02, 6.591e+01, 1.413e+01, 1.060e+00, 1.000e-02, 0.000e+00, 1.850e+00]))
        np.testing.assert_almost_equal(new_atomic_linelist['theoretical_depth'][:10], np.array([0.01 , 0.04 , 0.   , 0.   , 0.82 , 0.31 , 0.02 , 0.   , 0.959, 0.04 ]))



