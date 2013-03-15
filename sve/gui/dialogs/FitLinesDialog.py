import numpy as np
from CustomDialog import *

class FitLinesDialog(CustomDialog):
    def __init__(self, parent, title, vel_atomic, vel_telluric):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity respect to atomic lines (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_atomic
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity respect to telluric lines (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_telluric
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = ["VALD.300_1100nm"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel_telluric=None, updated_vel_atomic=None):
        self.results = None
        if updated_vel_atomic != None:
            self.__components[0]["default"] = updated_vel_atomic
        if updated_vel_telluric != None:
            self.__components[1]["default"] = updated_vel_telluric
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


