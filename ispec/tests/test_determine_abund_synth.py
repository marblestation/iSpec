import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestDetermineAbundSynth(unittest.TestCase):

    def test_determine_abundances_using_synth_spectra(self):
        code = "spectrum"
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
        #--- Model spectra ----------------------------------------------------------
        # Parameters
        initial_teff = 5771.0
        initial_logg = 4.44
        initial_MH = 0.00
        initial_alpha = 0.00
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = 1.60 # Sun
        initial_limb_darkening_coeff = 0.6
        initial_R = to_resolution
        initial_vrad = 0
        max_iterations = 6

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
        atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(star_spectrum['waveobs']), wave_top=np.max(star_spectrum['waveobs']))
        atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

        isotopes = ispec.read_isotope_data(isotope_file)



        # Load model atmospheres
        modeled_layers_pack = ispec.load_modeled_layers_pack(model)

        # Load SPECTRUM abundances
        solar_abundances = ispec.read_solar_abundances(solar_abundances_file)


        # Free parameters
        #free_params = ["teff", "logg", "MH", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
        free_params = ["vrad"]
        #free_params = []

        # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
        chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
        chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

        element_name = "Ca"
        free_abundances = ispec.create_free_abundances_structure([element_name], chemical_elements, solar_abundances)
        free_abundances['Abund'] += initial_MH # Scale to metallicity

        linelist_free_loggf = None

        # Line regions
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        # Select only the lines to get abundances from
        line_regions = line_regions[np.logical_or(line_regions['note'] == element_name+' 1', line_regions['note'] == element_name+' 2')]
        line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)

        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or we can create the segments on the fly:
        segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=True, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = False, \
                max_iterations=max_iterations, \
                tmp_dir = None, \
                code=code)

        expected_params = {
                            'MH': 0.0,
                            'R': 47000.0,
                            'alpha': 0.0,
                            'limb_darkening_coeff': 0.6,
                            'logg': 4.44,
                            'teff': 5771.0,
                            'vmac': 4.19,
                            'vmic': 1.07,
                            'vrad0000': -0.014236572551593383,
                            'vrad0001': -0.060747762730039546,
                            'vrad0002': -0.015112973054249907,
                            'vrad0003': -0.05478978560293956,
                            'vrad0004': -0.15687878317415702,
                            'vrad0005': -0.22541956406770822,
                            'vrad0006': -0.20350778494525035,
                            'vrad0007': -0.16016986490598595,
                            'vsini': 1.6
                           }
        for k, v in list(expected_params.items()):
            self.assertAlmostEqual(params[k], v)
        expected_errors = {
                            'MH': 0.0,
                            'R': 0.0,
                            'alpha': 0.0,
                            'limb_darkening_coeff': 0.0,
                            'logg': 0.0,
                            'teff': 0.0,
                            'vmac': 0.0,
                            'vmic': 0.0,
                            'vrad0000': 2.2869694859069014,
                            'vrad0001': 0.6752670286020491,
                            'vrad0002': 0.27650862905582735,
                            'vrad0003': 1.078220425624207,
                            'vrad0004': 0.47623752166713595,
                            'vrad0005': 0.2582966608244762,
                            'vrad0006': 0.46080252531829013,
                            'vrad0007': 0.23844625288619709,
                            'vsini': 0.0
                          }
        for k, v in list(expected_errors.items()):
            self.assertAlmostEqual(errors[k], v)
        self.assertEqual(len(stats_linemasks), 8)
        self.assertEqual(abundances_found['element'][0], 'Ca')
        self.assertAlmostEqual(abundances_found['[X/H]'][0], -0.003417007071945477)
        self.assertAlmostEqual(abundances_found['[X/Fe]'][0], -0.003417007071945477)
        self.assertAlmostEqual(abundances_found['e[X/H]'][0], 0.02708655339373642)
        self.assertAlmostEqual(abundances_found['e[X/Fe]'][0], 0.02708655339373642)
