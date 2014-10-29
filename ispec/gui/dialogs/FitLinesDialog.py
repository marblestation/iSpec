import numpy as np
from CustomDialog import *

class FitLinesDialog(CustomDialog):
    def __init__(self, parent, title, resolution, vel_telluric):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "Entry"
        component["text"] = "Resolution"
        component["text-type"] = "int" # float, int or str
        component["default"] = resolution
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
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = ["VALD_atom.300_1100nm", \
                "GESv4_atom_hfs_iso.475_685nm", "GESv4_atom_nohfs_noiso.475_685nm",  \
                "GESv4_atom_hfs_iso.845_895nm", "GESv4_atom_nohfs_noiso.845_895nm",  \
                "GESv5_atom_hfs_iso.420_920nm", "GESv5_atom_nohfs_noiso.420_920nm",  \
                "SEPv1.655_1020nm",  \
                "Kurucz_atom.300_1100nm", "NIST_atom.300_1100nm", "SPECTRUM.300_1000nm"]
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum atomic wavelength difference"
        component["text-type"] = "float" # float, int or str
        component["default"] = 0.005
        component["minvalue"] = 0.0001
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Allow peak position adjustment"
        component["default"] = True
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Check derivatives before fitting"
        component["default"] = False
        self.__components.append(component)

    def show(self, updated_vel_telluric=None):
        self.results = None
        if updated_vel_telluric is not None:
            self.__components[1]["default"] = updated_vel_telluric
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


