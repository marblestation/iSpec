import numpy as np
from CustomDialog import *

class EstimateSNRDialog(CustomDialog):
    def __init__(self, parent, title, num_points=10, resampling=0.001):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Label"
        component["text"] = "* Number of points and wavelength step\n is only used if SNR is estimated from fluxes (not errors)"
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Number of points"
        component["text-type"] = "int" # float, int or str
        component["default"] = num_points
        component["minvalue"] = 1.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength step (resampling)"
        component["text-type"] = "float" # float, int or str
        component["default"] = resampling
        component["minvalue"] = 1.0e-10
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Estimate SNR"
        component["options"] = ["Directly from reported errors", "From fluxes in blocks of N points"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


