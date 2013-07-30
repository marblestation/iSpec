import numpy as np
from CustomDialog import *

class CleanTelluricsDialog(CustomDialog):
    def __init__(self, parent, title, vel_telluric, min_vel, max_vel, min_depth):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity relative to tellurics"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_telluric
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Minimum velocity"
        component["text-type"] = "float" # float, int or str
        component["default"] = min_vel
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum velocity"
        component["text-type"] = "float" # float, int or str
        component["default"] = max_vel
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Minimum tellurics depth"
        component["text-type"] = "float" # float, int or str
        component["default"] = min_depth
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Replace by"
        #component["options"] = ["Zeros", "NaN", "Continuum", "Completely remove"]
        component["options"] = ["Zeros", "Continuum", "Completely remove"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel=None):
        self.results = None
        if updated_vel is not None:
            self.__components[0]["default"] = updated_vel
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


