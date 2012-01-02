import numpy as np
from common import *

# Radial vel in km/s
def correct_radial_velocity(spectra, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    # Correct wavelength scale for radial velocity
    spectra['waveobs'] = spectra['waveobs'] / ((radial_vel / c) + 1)
    return spectra

# Shift regions considering the known radial velocity of the star
def correct_radial_velocity_regions(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000
    # Oposite sense because we want to correct regions, not the spectra
    radial_vel = radial_vel * -1

    continuum_regions = read_continuum_regions(continuum_regions_filename)
    continuum_regions['wave_base'] = continuum_regions['wave_base'] / ((radial_vel / c) + 1)
    continuum_regions['wave_top'] = continuum_regions['wave_top'] / ((radial_vel / c) + 1)
    write_continuum_regions(continuum_regions, continuum_regions_filename_out)
     

    # Lines
    line_regions = read_line_regions(line_regions_filename)
    line_regions['wave_base'] = line_regions['wave_base'] / ((radial_vel / c) + 1)
    line_regions['wave_top'] = line_regions['wave_top'] / ((radial_vel / c) + 1)
    line_regions['wave_peak'] = line_regions['wave_peak'] / ((radial_vel / c) + 1)
    write_line_regions(line_regions, line_regions_filename_out)

    segment_regions = read_segment_regions(segment_regions_filename)
    segment_regions['wave_base'] = segment_regions['wave_base'] / ((radial_vel / c) + 1)
    segment_regions['wave_top'] = segment_regions['wave_top'] / ((radial_vel / c) + 1)
    write_segment_regions(segment_regions, segment_regions_filename_out)


#~ if __name__ == '__main__':
    #~ continuum_regions_filename = "input/LUMBA/UVES_MRD_sun_cmask.txt"
    #~ line_regions_filename = "input/LUMBA/UVES_MRD_sun_Fe-linelist.txt"
    #~ segment_regions_filename = "input/LUMBA/UVES_MRD_sun_segments.txt"
#~ 
    #~ radial_vel = -97.2 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+cas+a
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPD_mu_cas_a_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 14.03 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+leo
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 0 # The original spectra for Arcturus has been corrected for radial velocity
    #~ #radial_vel = -5.19 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=arcturus
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)











