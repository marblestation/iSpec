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
from .CustomDialog import *

class DegradeResolutionDialog(CustomDialog):
    def __init__(self, parent, title, from_resolution, to_resolution):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Label"
        component["text"] = "If initial resolution is set to zero, \nspectrum will be just smoothed by a gaussian \ndetermined by the final resolution \n(FWHM = wavelengths / final resolution)"
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Initial resolution"
        component["text-type"] = "float" # float, int or str
        component["default"] = from_resolution
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Final resolution"
        component["text-type"] = "float" # float, int or str
        component["default"] = to_resolution
        component["minvalue"] = 1.0
        component["maxvalue"] = np.inf
        self.__components.append(component)

    def show(self, updated_from_resolution=None, updated_to_resolution=None):
        self.results = None
        if updated_from_resolution is not None:
            self.__components[1]["default"] = updated_from_resolution
        if updated_to_resolution is not None:
            self.__components[2]["default"] = updated_to_resolution
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


