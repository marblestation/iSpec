import numpy as np
from CustomDialog import *

class SolverDialog(CustomDialog):
    def __init__(self, parent, title, resolution, teff, logg, feh, vmic, vmac, vsini, limb_darkening_coeff):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Model atmosphere"
        component["options"] = ["MARCS", "Kurucz", "Castelli"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = ["VALD.300_1100nm", "Kurucz.300_1100nm", "NIST.300_1100nm", "SPECTRUM.300_1000nm"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Effective temperature (K)"
        component["text-type"] = "float" # float, int or str
        component["default"] = teff
        component["minvalue"] = 2500
        component["maxvalue"] = 8000
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free Teff"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Surface gravity (log g)"
        component["text-type"] = "float" # float, int or str
        component["default"] = logg
        component["minvalue"] = 0
        component["maxvalue"] = 5
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free Log(g)"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Metallicity [Fe/H]"
        component["text-type"] = "float" # float, int or str
        component["default"] = feh
        component["minvalue"] = -5
        component["maxvalue"] = 1
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free [Fe/H]"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Microturbulence velocity (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vmic
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free Vmic"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Macroturbulence velocity (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vmac
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free Vmac"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Rotation (v sin(i)) (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vsini
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free vsin(i)"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Limb darkening coefficient"
        component["text-type"] = "float" # float, int or str
        component["default"] = limb_darkening_coeff
        component["minvalue"] = 0
        component["maxvalue"] = 1
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free limb dark. coeff."
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Resolution"
        component["text-type"] = "float" # float, int or str
        component["default"] = resolution
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free resolution"
        component["default"] = False
        self.__components.append(component)


    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


