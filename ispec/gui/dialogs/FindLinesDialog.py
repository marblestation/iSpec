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

class FindLinesDialog(CustomDialog):
    def __init__(self, parent, title, min_depth, max_depth, vel_telluric, resolution, elements, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Minimum depth (% of the continuum)"
        component["text-type"] = "float" # float, int or str
        component["default"] = min_depth
        component["minvalue"] = 0.0
        component["maxvalue"] = 1.0
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum depth (% of the continuum)"
        component["text-type"] = "float" # float, int or str
        component["default"] = max_depth
        component["minvalue"] = 1.0e-10
        component["maxvalue"] = 1.0
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Select elements (comma separated)"
        component["text-type"] = "str" # float, int or str
        component["default"] = elements
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Resolution"
        component["text-type"] = "int" # float, int or str
        component["default"] = resolution
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity respect to telluric lines (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_telluric
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Discard affected by tellurics"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Check derivatives before fitting"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = lists['atomic_lines']['name']
        component["default"] = component["options"][default_lists['atomic_lines']]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum atomic wavelength difference"
        component["text-type"] = "float" # float, int or str
        component["default"] = 0.005
        component["minvalue"] = 0.00001
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Look for line masks in"
        component["options"] = ["The whole spectra", "Only inside segments"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel_telluric=None):
        self.results = None
        if updated_vel_telluric is not None:
            self.__components[4]["default"] = updated_vel_telluric
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


