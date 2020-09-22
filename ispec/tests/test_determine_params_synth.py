import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestDetermineParamsSynth(unittest.TestCase):

    def test_determine_astrophysical_parameters_using_grid(self):
        star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Sun_Vesta-1.txt.gz")
        #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_Arcturus.txt.gz")
        #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muCas.txt.gz")
        #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/NARVAL_muLeo.txt.gz")
        #star_spectrum = ispec.read_spectrum(ispec_dir + "/input/spectra/examples/HARPS.GBOG_Procyon.txt.gz")
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
        initial_teff = 5750.0
        initial_logg = 4.5
        initial_MH = 0.00
        initial_alpha = 0.00
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = 2.0
        initial_limb_darkening_coeff = 0.6
        initial_R = to_resolution
        initial_vrad = 0
        max_iterations = 20

        code = "grid"
        precomputed_grid_dir = ispec_dir + "/input/grid/SPECTRUM_MARCS.GES_GESv6_atom_hfs_iso.480_680nm_light/"

        atomic_linelist = None
        isotopes = None
        modeled_layers_pack = None
        solar_abundances = None
        free_abundances = None
        linelist_free_loggf = None

        # Free parameters (vmic cannot be used as a free parameter when using a spectral grid)
        #free_params = ["teff", "logg", "MH", "alpha", "vmic", "vmac", "vsini", "R", "vrad", "limb_darkening_coeff"]
        free_params = ["teff", "logg", "MH", "alpha", "vmic", "R"]

        # Line regions
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        ## Select only some lines to speed up the execution (in a real analysis it is better not to do this)
        #line_regions = line_regions[np.logical_or(line_regions['note'] == 'Ti 1', line_regions['note'] == 'Ti 2')]
        #line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)
        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or we can create the segments on the fly:
        segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

        ## Add also regions from the wings of strong lines:
        # H beta
        hbeta_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Hbeta.txt")
        hbeta_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Hbeta_segments.txt")
        line_regions = np.hstack((line_regions, hbeta_lines))
        segments = np.hstack((segments, hbeta_segments))
        # H alpha
        halpha_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Halpha.txt")
        halpha_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Halpha_segments.txt")
        line_regions = np.hstack((line_regions, halpha_lines))
        segments = np.hstack((segments, halpha_segments))
        # Magnesium triplet
        mgtriplet_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_MgTriplet.txt")
        mgtriplet_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_MgTriplet_segments.txt")
        line_regions = np.hstack((line_regions, mgtriplet_lines))
        segments = np.hstack((segments, mgtriplet_segments))

        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=False, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = True, \
                max_iterations=max_iterations, \
                tmp_dir = None, \
                code=code, precomputed_grid_dir=precomputed_grid_dir)

        expected_params = {
                            'teff': 5851.628742661467,
                            'logg': 4.461276199140998,
                            'MH': 0.024255109120437632,
                            'alpha': 0.0,
                            'vmic': 0.5246819025206294,
                            'vmac': 4.51,
                            'vsini': 2.0,
                            'limb_darkening_coeff': 0.6,
                            'R': 77632.52714857543
        }
        for k, v in list(expected_params.items()):
            self.assertAlmostEqual(params[k], v)
        expected_errors = {
                            'teff': 14.225900341961772,
                            'logg': 0.02712239389209893,
                            'MH': 0.005795561724255147,
                            'alpha': 0.0,
                            'vmic': 0.042710575599214096,
                            'vmac': 0.0,
                            'vsini': 0.0,
                            'limb_darkening_coeff': 0.0,
                            'R': 3520.9717057141356
        }
        for k, v in list(expected_errors.items()):
            self.assertAlmostEqual(errors[k], v)
        self.assertEqual(len(stats_linemasks), 350)

    def test_determine_astrophysical_parameters_using_synth_spectra(self):
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
        initial_teff = 5750.0
        initial_logg = 4.5
        initial_MH = 0.00
        initial_alpha = ispec.determine_abundance_enchancements(initial_MH)
        initial_vmic = ispec.estimate_vmic(initial_teff, initial_logg, initial_MH)
        initial_vmac = ispec.estimate_vmac(initial_teff, initial_logg, initial_MH)
        initial_vsini = 2.0
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
        free_params = ["teff", "logg", "MH", "vmic", "R"]

        # Free individual element abundance
        free_abundances = None
        linelist_free_loggf = None

        # Line regions
        line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all.txt".format(code))
        #line_regions = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        ## Select only some lines to speed up the execution (in a real analysis it is better not to do this)
        line_regions = line_regions[np.logical_or(line_regions['note'] == 'Ti 1', line_regions['note'] == 'Ti 2')]
        line_regions = ispec.adjust_linemasks(normalized_star_spectrum, line_regions, max_margin=0.5)
        # Read segments if we have them or...
        #segments = ispec.read_segment_regions(ispec_dir + "/input/regions/fe_lines_segments.txt")
        # ... or we can create the segments on the fly:
        segments = ispec.create_segments_around_lines(line_regions, margin=0.25)

        ### Add also regions from the wings of strong lines:
        ## H beta
        #hbeta_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Hbeta.txt")
        #hbeta_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Hbeta_segments.txt")
        #line_regions = np.hstack((line_regions, hbeta_lines))
        #segments = np.hstack((segments, hbeta_segments))
        ## H alpha
        #halpha_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_Halpha.txt")
        #halpha_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_Halpha_segments.txt")
        #line_regions = np.hstack((line_regions, halpha_lines))
        #segments = np.hstack((segments, halpha_segments))
        ## Magnesium triplet
        #mgtriplet_lines = ispec.read_line_regions(ispec_dir + "input/regions/wings_MgTriplet.txt")
        #mgtriplet_segments = ispec.read_segment_regions(ispec_dir + "input/regions/wings_MgTriplet_segments.txt")
        #line_regions = np.hstack((line_regions, mgtriplet_lines))
        #segments = np.hstack((segments, mgtriplet_segments))

        obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                ispec.model_spectrum(normalized_star_spectrum, star_continuum_model, \
                modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                linemasks=line_regions, \
                enhance_abundances=True, \
                use_errors = True, \
                vmic_from_empirical_relation = False, \
                vmac_from_empirical_relation = True, \
                max_iterations=max_iterations, \
                tmp_dir = None, \
                code=code)

        expected_params = {
                            'teff': 5700.60895591792,
                            'logg': 4.357052646691409,
                            'MH': -0.10729296408531262,
                            'alpha': 0.042917185634125055,
                            'vmic': 1.1412534073088432,
                            'vmac': 4.04,
                            'vsini': 2.0,
                            'limb_darkening_coeff': 0.6,
                            'R': 50101.842351606625
        }
        for k, v in list(expected_params.items()):
            self.assertAlmostEqual(params[k], v)
        expected_errors = {
                            'teff': 66.2711106842635,
                            'logg': 0.10477704697086701,
                            'MH': 0.07756514816665908,
                            'alpha': 0.0,
                            'vmic': 0.11389528849649228,
                            'vmac': 0.0,
                            'vsini': 0.0,
                            'limb_darkening_coeff': 0.0,
                            'R': 3742.908365508009
        }
        for k, v in list(expected_errors.items()):
            self.assertAlmostEqual(errors[k], v)
        self.assertEqual(len(stats_linemasks), 29)
