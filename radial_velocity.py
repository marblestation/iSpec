"""
    This file is part of Spectra.
    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
    
    Spectra is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spectra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with Spectra.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
from common import *
from interpolate import *
from fitting import *
from pymodelfit import GaussianModel


## Determines the radial velocity profile using the lines regions:
## - For each point of each line region, get the spectra flux
## - The result of adding all these spectra fluxes will be a noisy base and a deep gaussian
## Return the radial velocity x coordenates and the averaged fluxes
def build_radial_velocity_profile(spectra, continuum, lines, rv_limit = 200, rv_step=0.5, frame=None):
    rv = 0 # km/s
    
    min_wave = np.min(spectra['waveobs'])
    max_wave = np.max(spectra['waveobs'])
    
    # Only use the lines that are on the spectra wavelength
    wfilter = (lines['wave_peak'] >= min_wave) & (lines['wave_peak'] <= max_wave)
    lines = lines[wfilter]
    
    if frame != None:
        frame.update_progress(0)
        total_lines = len(lines)
    
    # Speed of light in m/s
    c = 299792458.0
    
    # Build a window with the maximum size as to store the biggest region
    total_window = (2*rv_limit) / rv_step
    total_window = int(np.ceil(total_window))
    fluxes = np.zeros(total_window)
    
    line_num = 0
    # For each point of each line region, get the flux of the spectra and sum it all
    for line in lines:
        i = 0
        sindex = 0
        # Common sized window in km/s for all lines (but different in wavelength)
        window_base = line['wave_peak'] + (line['wave_peak'] * ((-1*rv_limit*1000) / c)) # nm
        window_top = line['wave_peak'] + (line['wave_peak'] * ((rv_limit*1000) / c))   # nm
        increment = (window_top - window_base) / total_window
        
        wavelength = window_base
        while wavelength <= window_top and i < total_window:
            # Outside spectra
            if wavelength > max_wave or wavelength < min_wave:
                fluxes[i] += 0 # Outside spectra
            else:
                flux, sindex = get_flux(spectra, wavelength)
                fluxes[i] += flux - continuum(wavelength) # Normalize continuum
            #print wavelength
            wavelength += increment
            i += 1
        if frame != None:
            current_work_progress = (line_num*1.0 / total_lines) * 100
            frame.update_progress(current_work_progress)
        line_num += 1
    
    fluxes = fluxes/line_num
    xcoord = np.linspace( -1*rv_limit, rv_limit, total_window ) # z numbers from x to y
    
    return xcoord, fluxes


## Fits a Gaussian to the radial velocity profile
## Return the fitted model
def model_radial_velocity_profile(xcoord, fluxes, renormalize=False):
    if renormalize:
        # Renormalize continuum
        renorm_model = UniformKnotSplineModel(nknots=1)
        renorm_model.fitData(xcoord, fluxes)
        fluxes = fluxes - renorm_model(xcoord)
    
    # The added flux should have a big gaussian
    model = GaussianModel()
    
    # First estimate for the gaussian mean: the place with the min flux
    model.mu = xcoord[np.argmin(fluxes)]
    model.sig = 20
    model.A = -0.70
    model.fitData(xcoord, fluxes, fixedpars=[])
    
    return model

## Constructs the radial velocity profile, fits a Gaussian and returns the mean (km/s)
def calculate_radial_velocity(spectra, continuum, lines, rv_limit = 200, rv_step=0.5, renormalize=False, frame=None):
    xcoord, fluxes = build_radial_velocity_profile(spectra, continuum, lines, rv_limit=rv_limit, rv_step=rv_step, frame=frame)
    model = model_radial_velocity_profile(xcoord, fluxes, renormalize=renormalize)

    rv = np.round(model.mu, 2)  # km/s
    ##rms = np.mean(model.residuals()) + np.std(model.residuals())
    ## Residual peak should be at least below 1 to be confident about the fit
    #residual_peak = model(model.mu) - np.min(fluxes)
    return rv

# Radial vel in km/s
def correct_radial_velocity(spectra, radial_vel):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    # Correct wavelength scale for radial velocity
    spectra['waveobs'] = spectra['waveobs'] / ((radial_vel / c) + 1)
    return spectra

def correct_radial_velocity_regions(regions, radial_vel, with_peak=False):
    # Speed of light in m/s
    c = 299792458.0
    # Radial velocity from km/s to m/s
    radial_vel = radial_vel * 1000

    regions['wave_base'] = regions['wave_base'] / ((radial_vel / c) + 1)
    regions['wave_top'] = regions['wave_top'] / ((radial_vel / c) + 1)
    if with_peak:
        regions['wave_peak'] = regions['wave_peak'] / ((radial_vel / c) + 1)
    return regions

# Shift regions considering the known radial velocity of the star
def correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel):
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
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 14.03 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=mu+leo
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MRG_mu_leo_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)
#~ 
    #~ radial_vel = 0 # The original spectra for Arcturus has been corrected for radial velocity
    #~ #radial_vel = -5.19 # km/s http://simbad.u-strasbg.fr/simbad/sim-basic?Ident=arcturus
#~ 
    #~ continuum_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_cmask.txt"
    #~ line_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_Fe-linelist.txt"
    #~ segment_regions_filename_out = "input/LUMBA/UVES_MPG_arcturus_segments.txt"
    #~ 
    #~ correct_radial_velocity_regions_files(continuum_regions_filename, line_regions_filename, segment_regions_filename, continuum_regions_filename_out, line_regions_filename_out, segment_regions_filename_out, radial_vel)











