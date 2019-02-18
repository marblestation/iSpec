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

class FitContinuumDialog(CustomDialog):
    def __init__(self, parent, title, R, nknots, degrees, median_wave_range, max_wave_range, strong_line_probability, templates, fit_type):
        self.__parent = parent
        self.__title = title
        self.__components = []
        self.__fit_type = fit_type

        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Fitting model"
        #component["options"] = ["Splines", "Polynomy", "Template", "Fixed value"]
        component["options"] = [fit_type]
        component["default"] = component["options"][0]
        self.__components.append(component)

        if fit_type == "Splines":
            component = {}
            component["type"] = "Label"
            component["text"] = "Suggested number of splines: %i" % nknots
            self.__components.append(component)
            component = {}
            component["type"] = "Entry"
            component["text"] = "Number of splines"
            component["text-type"] = "int" # float, int or str
            component["default"] = nknots
            component["minvalue"] = 1.0
            component["maxvalue"] = np.inf
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy']:
            component = {}
            component["type"] = "Entry"
            component["text"] = "Degree"
            component["text-type"] = "int" # float, int or str
            component["default"] = degrees
            component["minvalue"] = 1.0
            component["maxvalue"] = np.inf
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy', 'Template']:
            component = {}
            component["type"] = "Entry"
            component["text"] = "Resolution"
            component["text-type"] = "int" # float, int or str
            component["default"] = R
            component["minvalue"] = 0.0
            component["maxvalue"] = np.inf
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy']:
            component = {}
            component["type"] = "OptionMenu"
            component["text"] = "Filtering order"
            component["options"] = ["median+max", "max+median"]
            component["default"] = component["options"][0]
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy', 'Template']:
            component = {}
            component["type"] = "Entry"
            component["text"] = "Wavelength step for median selection"
            component["text-type"] = "float" # float, int or str
            component["default"] = median_wave_range
            component["minvalue"] = 0.0
            component["maxvalue"] = np.inf
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy']:
            component = {}
            component["type"] = "Entry"
            component["text"] = "Wavelength step for max selection"
            component["text-type"] = "float" # float, int or str
            component["default"] = max_wave_range
            component["minvalue"] = 0.0
            component["maxvalue"] = np.inf
            self.__components.append(component)
            component = {}
            component["type"] = "Checkbutton"
            component["text"] = "Use spectrum's errors as weights for the fitting process"
            component["default"] = True
            self.__components.append(component)
            component = {}
            component["type"] = "Checkbutton"
            component["text"] = "Automatically find and ignore strong lines"
            component["default"] = True
            self.__components.append(component)
            component = {}
            component["type"] = "Entry"
            component["text"] = "Strong line probability threshold"
            component["text-type"] = "float" # float, int or str
            component["default"] = strong_line_probability
            component["minvalue"] = 0.0
            component["maxvalue"] = 1.0
            self.__components.append(component)

        if fit_type in ['Splines', 'Polynomy', 'Template']:
            component = {}
            component["type"] = "Checkbutton"
            component["text"] = "Consider only continuum regions"
            component["default"] = False
            self.__components.append(component)
            component = {}
            component["type"] = "Checkbutton"
            component["text"] = "Ignore line regions"
            component["default"] = False
            self.__components.append(component)
            component = {}
            component["type"] = "Checkbutton"
            component["text"] = "Treat each segment independently"
            component["default"] = False
            self.__components.append(component)

        if fit_type in ['Template']:
            component = {}
            component["type"] = "OptionMenu"
            component["text"] = "Use as a template"
            component["options"] = templates
            component["default"] = templates[0]
            self.__components.append(component)

        if fit_type in ['Fixed value']:
            component = {}
            component["type"] = "Entry"
            component["text"] = "Fixed value"
            component["text-type"] = "float" # float, int or str
            component["default"] = 1.0
            component["minvalue"] = -np.inf
            component["maxvalue"] = np.inf
            self.__components.append(component)

    def show(self, suggested_nknots=None, updated_templates=None):
        self.results = None
        if self.__fit_type in ['Splines'] and suggested_nknots is not None:
            self.__components[1]["text"] = "Suggested number of splines\n based on the wavelength range: %i" % suggested_nknots

        if self.__fit_type == "Template" and updated_templates is not None:
            self.__components[6]["options"] = updated_templates
            # Validate that the default value (previous user selected value) exists in the new template list
            default_ok = False
            for template in updated_templates:
                if template == self.__components[6]["default"]:
                    default_ok = True
            if not default_ok:
                self.__components[6]["default"] = updated_templates[0]
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


