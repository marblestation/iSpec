import numpy as np
from CustomDialog import *

class FitContinuumDialog(CustomDialog):
    def __init__(self, parent, title, nknots, median_wave_range, max_wave_range):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Label"
        component["text"] = "Suggested number of polynomial degrees/splines: %i" % nknots
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Number of polynomial degrees/splines"
        component["text-type"] = "int" # float, int or str
        component["default"] = nknots
        component["minvalue"] = 1.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength step for median selection"
        component["text-type"] = "float" # float, int or str
        component["default"] = median_wave_range
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength step for max selection"
        component["text-type"] = "float" # float, int or str
        component["default"] = max_wave_range
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Fit using"
        component["options"] = ["The whole spectra", "Only continuum regions"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Fitting model"
        component["options"] = ["Polynomy", "Splines", "Fixed value"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Fixed value"
        component["text-type"] = "float" # float, int or str
        component["default"] = 1.0
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)

    def show(self, suggested_nknots=None):
        self.results = None
        if suggested_nknots != None:
            self.__components[0]["text"] = "Suggested number of polynomial degrees/splines\n based on the wavelength range: %i" % suggested_nknots
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


