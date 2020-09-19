import os
import sys
import unittest
import numpy as np

ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


class TestDetermineLoggfSynth(unittest.TestCase):

    def test_determine_loggf_line_by_line_using_synth_spectra(self):
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
        #free_params = ["vrad"]
        free_params = []

        # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
        chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
        chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

        # Line regions
        line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_GES/{}_synth_good_for_params_all_extended.txt".format(code))
        #line_regions_with_atomic_data = ispec.read_line_regions(ispec_dir + "/input/regions/47000_VALD/{}_synth_good_for_params_all_extended.txt".format(code))
        # Select only the lines to get abundances from
        line_regions_with_atomic_data = line_regions_with_atomic_data[:5]
        line_regions_with_atomic_data = ispec.adjust_linemasks(normalized_star_spectrum, line_regions_with_atomic_data, max_margin=0.5)

        output_dirname = "example_loggf_line_by_line_%s" % (code,)
        #ispec.mkdir_p(output_dirname)
        for i, line in enumerate(line_regions_with_atomic_data):
            # Directory and file names
            #element_name = "_".join(line['element'].split())
            element_name = "_".join(line['note'].split())
            common_filename = "example_" + code + "_individual_" + element_name + "_%.4f" % line['wave_peak']

            # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
            free_abundances = None

            # Line by line
            individual_line_regions = line_regions_with_atomic_data[i:i+1] # Keep recarray structure
            linelist_free_loggf = line_regions_with_atomic_data[i:i+1] # Keep recarray structure

            # Filter the line that we want to determine the loggf from the global atomic linelist
            lfilter = atomic_linelist['element'] == linelist_free_loggf['element'][0]
            for key in ['wave_nm', 'lower_state_eV', 'loggf', 'stark', 'rad', 'waals']:
                lfilter = np.logical_and(lfilter, np.abs(atomic_linelist[key] - linelist_free_loggf[key][0]) < 1e-9)

            # Segment
            segments = ispec.create_segments_around_lines(individual_line_regions, margin=0.25)
            wfilter = ispec.create_wavelength_filter(normalized_star_spectrum, regions=segments) # Only use the segment

            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
                    ispec.model_spectrum(normalized_star_spectrum[wfilter], star_continuum_model, \
                    modeled_layers_pack, atomic_linelist[~lfilter], isotopes, solar_abundances, free_abundances, linelist_free_loggf, initial_teff, \
                    initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, \
                    initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=segments, \
                    linemasks=individual_line_regions, \
                    enhance_abundances=True, \
                    use_errors = True, \
                    vmic_from_empirical_relation = False, \
                    vmac_from_empirical_relation = False, \
                    max_iterations=max_iterations, \
                    tmp_dir = None, \
                    code=code)

            self.assertAlmostEqual(loggf_found['loggf'][0], -1.0794949759801327)
            self.assertAlmostEqual(loggf_found['eloggf'][0], 0.11403936780212467)
            self.assertEqual(len(loggf_found['loggf']), 1)
            self.assertEqual(individual_line_regions['lower_state_eV'][0], loggf_found['linelist']['lower_state_eV'][0])
            break

