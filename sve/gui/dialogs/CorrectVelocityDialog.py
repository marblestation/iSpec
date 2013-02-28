import numpy as np
from CustomDialog import *

class CorrectVelocityDialog(CustomDialog):
    def __init__(self, parent, title, vel_type, rv):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity relative to %s (km/s)" % vel_type
        component["text-type"] = "float" # float, int or str
        component["default"] = rv
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Apply correction on"
        component["options"] = ["Spectra", "Regions"]
        component["default"] = component["options"][0]
        self.__components.append(component)

    def show(self, updated_vel=None):
        self.results = None
        if updated_vel != None:
            self.__components[0]["default"] = updated_vel
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


