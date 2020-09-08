import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestSynthesis(unittest.TestCase):

    def test_synthesize_spectrum_with_spectrum(self):
        synth_spectrum = self._synthesize_spectrum(code="spectrum")
        np.testing.assert_almost_equal(synth_spectrum['flux'][:10], np.array([0.        , 0.99776054, 0.99698647, 0.99558204, 0.99346707, 0.99032837, 0.98571305, 0.9790391 , 0.96960094, 0.95660359]))
        np.testing.assert_almost_equal(synth_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))

    def test_synthesize_spectrum_with_turbospectrum(self):
        synth_spectrum = self._synthesize_spectrum(code="turbospectrum")
        np.testing.assert_almost_equal(synth_spectrum['flux'][:10], np.array([1.00000000e-10, 9.65410687e-01, 9.65625883e-01, 9.67300192e-01, 9.69917913e-01, 9.71878922e-01, 9.71699328e-01, 9.68304584e-01, 9.60961187e-01, 9.49084381e-01]))
        np.testing.assert_almost_equal(synth_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))

    def test_synthesize_spectrum_with_sme(self):
        synth_spectrum = self._synthesize_spectrum(code="sme")
        np.testing.assert_almost_equal(synth_spectrum['flux'][:10], np.array([1.00000000e-10, 9.66205259e-01, 9.66583367e-01, 9.68499306e-01, 9.71398442e-01, 9.73687695e-01, 9.73901693e-01, 9.70988650e-01, 9.64237832e-01, 9.53080035e-01]))
        np.testing.assert_almost_equal(synth_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))

    def test_synthesize_spectrum_with_moog(self):
        synth_spectrum = self._synthesize_spectrum(code="moog")
        np.testing.assert_almost_equal(synth_spectrum['flux'][:10], np.array([1.00000000e-10, 9.66309632e-01, 9.66671555e-01, 9.68536430e-01, 9.71328504e-01, 9.73431386e-01, 9.73350797e-01, 9.70000149e-01, 9.62629763e-01, 9.50632869e-01]))
        np.testing.assert_almost_equal(synth_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))

    def test_synthesize_spectrum_with_synthe(self):
        synth_spectrum = self._synthesize_spectrum(code="synthe")
        np.testing.assert_almost_equal(synth_spectrum['flux'][:10], np.array([1.00000000e-10, 9.66674590e-01, 9.67052538e-01, 9.68907744e-01, 9.71640919e-01, 9.73635560e-01, 9.73411213e-01, 9.69907277e-01, 9.62408637e-01, 9.50346662e-01]))
        np.testing.assert_almost_equal(synth_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))

    def _synthesize_spectrum(self, code):
        #--- Synthesizing spectrum -----------------------------------------------------
        # Parameters
        teff = 5771.0
        logg = 4.44
        MH = 0.00
        alpha = ispec.determine_abundance_enchancements(MH)
        microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
        macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
        vsini = 1.60 # Sun
        limb_darkening_coeff = 0.6
        resolution = 300000
        wave_step = 0.001

        # Wavelengths to synthesis
        #regions = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        regions = None
        wave_base = 515.0 # Magnesium triplet region
        wave_top = 525.0


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

        isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

        # Load chemical information and linelist
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

        isotopes = ispec.read_isotope_data(isotope_file)

        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
        #solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)
        # Load SPECTRUM abundances
        solar_abundances = ispec.read_solar_abundances(solar_abundances_file)

        ## Custom fixed abundances
        #fixed_abundances = ispec.create_free_abundances_structure(["C", "N", "O"], chemical_elements, solar_abundances)
        #fixed_abundances['Abund'] = [-3.49, -3.71, -3.54] # Abundances in SPECTRUM scale (i.e., x - 12.0 - 0.036) and in the same order ["C", "N", "O"]
        ## No fixed abundances
        fixed_abundances = None

        # Validate parameters
        if not ispec.valid_atmosphere_target(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                    fall out of theatmospheric models."
            print(msg)

        # Prepare atmosphere model
        atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)

        # Synthesis
        synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
        synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
                atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances, \
                fixed_abundances, microturbulence_vel = microturbulence_vel, \
                macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
                R=resolution, regions=regions, verbose=1,
                code=code)
        return synth_spectrum

    def test_interpolate_spectrum(self):
        #--- Synthesizing spectrum -----------------------------------------------------
        # Parameters
        #teff = 5771.0
        teff = 4025.0
        logg = 4.44
        MH = 0.00
        alpha = ispec.determine_abundance_enchancements(MH)
        microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) # 1.07
        macroturbulence = ispec.estimate_vmac(teff, logg, MH) # 4.21
        vsini = 1.60 # Sun
        limb_darkening_coeff = 0.6
        resolution = 300000
        wave_step = 0.001

        # Wavelengths to synthesis
        wave_base = 515.0 # Magnesium triplet region
        wave_top = 525.0

        code = "grid"
        precomputed_grid_dir = ispec_dir + "/input/grid/SPECTRUM_MARCS.GES_GESv5_atom_hfs_iso.480_680nm_light/"
        grid = ispec.load_spectral_grid(precomputed_grid_dir)

        atomic_linelist = None
        isotopes = None
        modeled_layers_pack = None
        solar_abundances = None
        fixed_abundances = None
        abundances = None
        atmosphere_layers = None
        regions = None

        # Validate parameters
        if not ispec.valid_interpolated_spectrum_target(grid, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic': microturbulence_vel}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                    fall out of the spectral grid limits."
            print(msg)

        # Interpolation
        interpolated_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
        interpolated_spectrum['flux'] = ispec.generate_spectrum(interpolated_spectrum['waveobs'], \
                atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, abundances, \
                fixed_abundances, microturbulence_vel = microturbulence_vel, \
                macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
                R=resolution, regions=regions, verbose=1,
                code=code, grid=grid)
        np.testing.assert_almost_equal(interpolated_spectrum['flux'][:10], np.array([1.00000000e-10, 9.29308407e-01, 9.21892388e-01, 9.11929100e-01, 9.01350821e-01, 8.89488008e-01, 8.75318514e-01, 8.57946158e-01, 8.36730900e-01, 8.11337422e-01]))
        np.testing.assert_almost_equal(interpolated_spectrum['waveobs'][:10], np.array([515.   , 515.001, 515.002, 515.003, 515.004, 515.005, 515.006, 515.007, 515.008, 515.009]))


