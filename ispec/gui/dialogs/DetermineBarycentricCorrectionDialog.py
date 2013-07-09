import numpy as np
from CustomDialog import *

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


