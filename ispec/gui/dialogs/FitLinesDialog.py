from __future__ import absolute_import
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

class FitLinesDialog(CustomDialog):
    def __init__(self, parent, title, resolution, vel_telluric, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__components = []
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
        component["minvalue"] = 0.0001
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Allow peak position adjustment"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Check derivatives before fitting"
        component["default"] = False
        self.__components.append(component)

    def show(self, updated_vel_telluric=None):
        self.results = None
        if updated_vel_telluric is not None:
            self.__components[1]["default"] = updated_vel_telluric
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


