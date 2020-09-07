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

class DetermineBarycentricCorrectionDialog(CustomDialog):
    def __init__(self, parent, title, date_string, time_string, ra, dec):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Date (DD/MM/YYY)"
        component["text-type"] = "str" # float, int or str
        component["default"] = date_string
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Time (HH:MM:SS)"
        component["text-type"] = "str" # float, int or str
        component["default"] = time_string
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Label"
        component["text"] = "Epoch J2000.0"
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Right ascension (HH:MM:SS)"
        component["text-type"] = "str" # float, int or str
        component["default"] = ra
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Declination (DD:MM:SS)"
        component["text-type"] = "str" # float, int or str
        component["default"] = dec
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


