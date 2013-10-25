import numpy as np
from CustomDialog import *

class SyntheticSpectrumDialog(CustomDialog):
    def __init__(self, parent, title, wave_base, wave_top, wave_step, resolution, teff, logg, feh, vmic, vmac, vsini, limb_darkening_coeff):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Model atmosphere"
        component["options"] = ["MARCS", "MARCS.GES", "MARCS.APOGEE", "ATLAS9.APOGEE", "ATLAS9.Castelli", "ATLAS9.Kurucz", "ATLAS9.Kirby"]
        component["default"] = component["options"][1]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Solar abundances"
        component["options"] = ["Asplund.2009", "Asplund.2005", "Grevesse.2007", "Grevesse.1998", "Anders.1989"]
        component["default"] = component["options"][2]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = ["VALD_atom.300_1100nm", \
                "GESv3.475_685nm", "GESv3_noABO.475_685nm", "GESv3_atom.475_685nm", "GESv3_atom_noABO.475_685nm", \
                "GESv4_atom.475_685nm", "GESv4_atom_noABO.475_685nm", "GESv4_atom_hfs.475_685nm", "GESv4_atom_hfs_noABO.475_685nm", \
                "GESv4_atom.845_895nm", "GESv4_atom_noABO.845_895nm", "GESv4_atom_hfs.845_895nm", "GESv4_atom_hfs_noABO.845_895nm", \
                "SEPv1.655_1020nm", "SEPv1_noABO.655_1020nm", \
                "Kurucz_atom.300_1100nm", "NIST_atom.300_1100nm", "SPECTRUM.300_1000nm"]
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
        component["type"] = "Entry"
        component["text"] = "Surface gravity (log g)"
        component["text-type"] = "float" # float, int or str
        component["default"] = logg
        component["minvalue"] = 0
        component["maxvalue"] = 5
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
        component["type"] = "Entry"
        component["text"] = "Microturbulence velocity (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vmic
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
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
        component["type"] = "Entry"
        component["text"] = "Rotation (v sin(i)) (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = vsini
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
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
        component["type"] = "Entry"
        component["text"] = "Resolution"
        component["text-type"] = "float" # float, int or str
        component["default"] = resolution
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength min (nm)"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_base
        component["minvalue"] = 300
        component["maxvalue"] = 1100
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength max (nm)"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_top
        component["minvalue"] = 300
        component["maxvalue"] = 1100
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Wavelength step (nm)"
        component["text-type"] = "float" # float, int or str
        component["default"] = wave_step
        component["minvalue"] = 1.0e-5
        component["maxvalue"] = 0.01
        self.__components.append(component)
        component = {}
        component["type"] = "Radiobutton"
        component["text"] = "Generate spectrum for"
        component["options"] = ["Custom range (defined above)", "Segments", "Line masks"]
        component["default"] = component["options"][0]
        self.__components.append(component)


    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


