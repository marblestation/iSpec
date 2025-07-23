#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
from ispec.abundances import enhance_solar_abundances
from . import moog
from . import sme
from . import spectrum
from . import synthe
from . import turbospectrum
from . import grid as grid_module

def generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, molecules_files=None, isotope_file=None, regions=None, code="spectrum", use_molecules=False, grid=None, tmp_dir=None):
    """
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'synthe', 'sme', 'grid']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    if code != "grid":
        # Filter out lines not supported by a given synthesizer
        if code == 'moog-scat':
            # MOOG-SCAT is backward compatible with MOOG
            linelist_code = 'moog'
        else:
            linelist_code = code
        lcode = linelist[linelist_code+'_support'] == "T"
        linelist = linelist[lcode]
        # Limit linelist to the region asked to be synthesized
        # Provide some margin or near-by deep lines might be omitted
        margin = 2. # 2 nm
        wfilter = np.logical_and(linelist['wave_nm'] >= np.min(waveobs)-margin, linelist['wave_nm'] <= np.max(waveobs)+margin)
        linelist = linelist[wfilter]

        abundances = enhance_solar_abundances(abundances, alpha)

        if code == "turbospectrum":
            return turbospectrum.generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)
        elif code in ("moog", "moog-scat"):
            return moog.generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, tmp_dir=tmp_dir, timeout=timeout, code=code)
        elif code == "synthe":
            return synthe.generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
        elif code == "sme":
            return sme.generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
        elif code == "spectrum":
            return spectrum.generate_fundamental_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions, tmp_dir=tmp_dir)
    elif code == "grid":
        if grid is None:
            raise Exception("Grid needed to generate an interpolated spectrum from a grid.")
        return grid_module.generate_fundamental_spectrum(grid, waveobs, teff, logg, MH, alpha, microturbulence_vel, regions=regions)


def generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = 2.0, macroturbulence = 3.0, vsini = 2.0, limb_darkening_coeff = 0.60, R=500000, verbose=0, gui_queue=None, timeout=1800, atmosphere_layers_file=None, abundances_file=None, fixed_abundances_file=None, linelist_file=None, molecules_files=None, isotope_file=None, regions=None, code="spectrum", use_molecules=False, grid=None, tmp_dir=None):
    """
    """
    code = code.lower()
    if code not in ['spectrum', 'turbospectrum', 'moog', 'moog-scat', 'synthe', 'sme', 'grid']:
        raise Exception("Unknown radiative transfer code: %s" % (code))

    if code != "grid":
        # Filter out lines not supported by a given synthesizer
        if code == 'moog-scat':
            # MOOG-SCAT is backward compatible with MOOG
            linelist_code = 'moog'
        else:
            linelist_code = code
        lcode = linelist[linelist_code+'_support'] == "T"
        linelist = linelist[lcode]
        # Limit linelist to the region asked to be synthesized
        # Provide some margin or near-by deep lines might be omitted
        margin = 2. # 2 nm
        wfilter = np.logical_and(linelist['wave_nm'] >= np.min(waveobs)-margin, linelist['wave_nm'] <= np.max(waveobs)+margin)
        linelist = linelist[wfilter]

        abundances = enhance_solar_abundances(abundances, alpha)

        if code == "turbospectrum":
            return turbospectrum.generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, regions=regions, use_molecules=use_molecules, tmp_dir=tmp_dir, timeout=timeout)
        elif code in ("moog", "moog-scat"):
            return moog.generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, regions=regions, tmp_dir=tmp_dir, timeout=timeout, code=code)
        elif code == "synthe":
            return synthe.generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose,  atmosphere_layers_file=atmosphere_layers_file, linelist_file=linelist_file, molecules_files=molecules_files, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
        elif code == "sme":
            return sme.generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose, regions=regions, tmp_dir=tmp_dir, timeout=timeout)
        elif code == "spectrum":
            return spectrum.generate_spectrum(waveobs, atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, verbose=verbose, gui_queue=gui_queue, timeout=timeout, atmosphere_layers_file=atmosphere_layers_file, abundances_file=abundances_file, fixed_abundances_file=fixed_abundances_file, linelist_file=linelist_file, isotope_file=isotope_file, regions=regions, tmp_dir=tmp_dir)
    elif code == "grid":
        if grid is None:
            raise Exception("Grid needed to generate an interpolated spectrum from a grid.")
        return grid_module.generate_spectrum(grid, waveobs, teff, logg, MH, alpha, microturbulence_vel, macroturbulence=macroturbulence, R=R, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, regions=regions)



