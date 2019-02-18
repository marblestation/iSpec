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
from CustomDialog import *

class CorrectVelocityDialog(CustomDialog):
    def __init__(self, parent, title, vel_type, rv):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity relative to %s (km/s)" % vel_type
        component["text-type"] = "float" # float, int or str
        component["default"] = rv
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Apply correction on"
        component["options"] = ["Spectra", "Regions"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel=None):
        self.results = None
        if updated_vel is not None:
            self.__components[0]["default"] = updated_vel
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


