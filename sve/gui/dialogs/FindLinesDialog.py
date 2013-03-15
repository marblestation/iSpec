import numpy as np
from CustomDialog import *

class FindLinesDialog(CustomDialog):
    def __init__(self, parent, title, min_depth=0.05, max_depth=1.0, vel_atomic=0.0, vel_telluric=0.0, resolution=300000, elements="Fe 1, Fe 2"):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Minimum depth (% of the continuum)"
        component["text-type"] = "float" # float, int or str
        component["default"] = min_depth
        component["minvalue"] = 0.0
        component["maxvalue"] = 1.0
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum depth (% of the continuum)"
        component["text-type"] = "float" # float, int or str
        component["default"] = max_depth
        component["minvalue"] = 1.0e-10
        component["maxvalue"] = 1.0
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Select elements (comma separated)"
        component["text-type"] = "str" # float, int or str
        component["default"] = elements
        component["minvalue"] = None
        component["maxvalue"] = None
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Resolution"
        component["text-type"] = "float" # float, int or str
        component["default"] = resolution
        component["minvalue"] = 1.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity respect to atomic lines (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_atomic
        component["minvalue"] = 0.0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity respect to telluric lines (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vel_telluric
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Discard affected by tellurics"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = ["VALD.300_1100nm"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Look for line masks in"
        component["options"] = ["The whole spectra", "Only inside segments"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel_atomic=None, updated_vel_telluric=None):
        self.results = None
        if updated_vel_atomic != None:
            self.__components[4]["default"] = updated_vel_atomic
        if updated_vel_telluric != None:
            self.__components[5]["default"] = updated_vel_telluric
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


