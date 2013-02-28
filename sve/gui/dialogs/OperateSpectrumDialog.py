import numpy as np
from CustomDialog import *

class OperateSpectrumDialog(CustomDialog):
    def __init__(self, parent, title, operations_description, waveobs, flux, err):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Available functions"
        component["options"] = operations_description
        component["default"] = component["options"][0]
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelengths ="
        component["text-type"] = "str" # float, int or str
        component["default"] = waveobs
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Fluxes ="
        component["text-type"] = "str" # float, int or str
        component["default"] = flux
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Errors ="
        component["text-type"] = "str" # float, int or str
        component["default"] = err
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


