import numpy as np
from CustomDialog import *

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
        if updated_from_resolution != None:
            self.__components[1]["default"] = updated_from_resolution
        if updated_to_resolution != None:
            self.__components[2]["default"] = updated_to_resolution
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


