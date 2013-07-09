import numpy as np
from CustomDialog import *

class CombineSpectraDialog(CustomDialog):
    def __init__(self, parent, title, wave_base, wave_top, wave_step):
        self.__parent = parent
        self.__title = title
        self.__components = []
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
        component["type"] = "Radiobutton"
        component["text"] = "Operation"
        component["options"] = ["Median", "Mean", "Subtract", "Add", "Divide"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


