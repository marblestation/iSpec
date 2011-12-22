import numpy as np

# Radial vel in km/s
def correct_radial_velocity(spectra, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    # Correct wavelength scale for radial velocity
    spectra['waveobs'] = spectra['waveobs'] / ((radial_vel / c) + 1)
    return spectra





