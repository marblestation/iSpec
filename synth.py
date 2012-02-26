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
import synthesizer
from atmospheres import *

modeled_layers_pack = load_modeled_layers_pack(filename='input/atmospheres/default.modeled_layers_pack.dump')

teff_obj = 5750.0
logg_obj = 4.5
MH_obj = 0.00

valid_objective(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
layers = interpolate_atmosphere_layers(modeled_layers_pack, teff_obj, logg_obj, MH_obj)
atm_filename = write_atmosphere(teff_obj, logg_obj, MH_obj, layers)


waveobs = np.arange(515.0, 520.0, 0.05)
fluxes = synthesizer.spectrum(waveobs*10.0, atm_filename, verbose=0)
print fluxes
os.remove(atm_filename)
