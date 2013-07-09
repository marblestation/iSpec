import numpy as np
from CustomDialog import *

class CleanSpectrumDialog(CustomDialog):
    def __init__(self, parent, title, flux_base, flux_top, err_base, err_top):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Filter by flux"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Base flux"
        component["text-type"] = "float" # float, int or str
        component["default"] = flux_base
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Top flux"
        component["text-type"] = "float" # float, int or str
        component["default"] = flux_top
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Filter by error"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Base error"
        component["text-type"] = "float" # float, int or str
        component["default"] = err_base
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Top error"
        component["text-type"] = "float" # float, int or str
        component["default"] = err_top
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


