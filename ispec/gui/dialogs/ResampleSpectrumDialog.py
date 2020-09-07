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

class ResampleSpectrumDialog(CustomDialog):
    def __init__(self, parent, title, wave_base, wave_top, wave_step, median_step, mean_step, min_step, max_step):
        self.__parent = parent
        self.__title = title
        self.__components = []

        self.__mean_step = mean_step
        self.__median_step = median_step
        self.__min_step = min_step
        self.__max_step = max_step
        component = {}
        component["type"] = "Entry"
        component["text"] = "Base wavelength"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_base
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Top wavelength"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_top
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength step"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_step
        component["minvalue"] = 1.0e-10
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Label"
        component["text"] = "Current wavelength step:\n min:" + str(self.__min_step) + " | max:" + str(self.__max_step) + "\nmedian:" + str(self.__median_step) + " | mean:" + str(self.__mean_step)
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Method"
        component["options"] = ["Linear", "Bessel", "Spline"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_mean_step=None, updated_median_step=None, updated_min_step=None, updated_max_step=None):
        self.results = None
        if updated_mean_step is not None:
            self.__mean_step = updated_mean_step
        if updated_median_step is not None:
            self.__median_step = updated_median_step
        if updated_min_step is not None:
            self.__min_step = updated_min_step
        if updated_max_step is not None:
            self.__max_step = updated_max_step
        self.__components[3]["text"] = "Current wavelength step:\n min:" + str(self.__min_step) + " | max:" + str(self.__max_step) + "\nmedian:" + str(self.__median_step) + " | mean:" + str(self.__mean_step)
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


