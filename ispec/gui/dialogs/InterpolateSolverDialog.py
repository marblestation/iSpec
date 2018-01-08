import numpy as np
from CustomDialog import *

class InterpolateSolverDialog(CustomDialog):
    def __init__(self, parent, title, resolution, teff, logg, feh, alpha, vmic, vmac, vsini, limb_darkening_coeff, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Grid"
        component["options"] = lists['grid']['name']
        component["default"] = component["options"][default_lists['grid']]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Effective temperature (K)"
        component["text-type"] = "float" # float, int or str
        component["default"] = teff
        component["minvalue"] = 2500
        component["maxvalue"] = 50000
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
        component["text"] = "Metallicity [M/H]"
        component["text-type"] = "float" # float, int or str
        component["default"] = feh
        component["minvalue"] = -5
        component["maxvalue"] = 1
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free [M/H]"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Alpha enhancement [alpha/Fe]"
        component["text-type"] = "float" # float, int or str
        component["default"] = alpha
        component["minvalue"] = -2
        component["maxvalue"] = 2
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free [alpha/Fe]"
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
        component["default"] = True
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
        component["type"] = "Checkbutton"
        component["text"] = "Automatic Vmac from empirical relation"
        component["default"] = True
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
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Radial velocity"
        component["text-type"] = "float" # float, int or str
        component["default"] = 0.
        component["minvalue"] = -5
        component["maxvalue"] = 5
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free radial velocity"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum number of iterations"
        component["text-type"] = "int" # float, int or str
        component["default"] = "6"
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)


    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


