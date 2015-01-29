import numpy as np
from CustomDialog import *

class SolverDialog(CustomDialog):
    def __init__(self, parent, title, resolution, teff, logg, feh, vmic, vmac, vsini, limb_darkening_coeff, lists, defaults_lists):
        self.__parent = parent
        self.__title = title
        self.__components = []
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Model atmosphere"
        component["options"] = lists['atmospheres']['name']
        component["default"] = component["options"][defaults_lists['atmospheres']]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Solar abundances"
        component["options"] = lists['abundances']['name']
        component["default"] = component["options"][defaults_lists['abundances']]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Line list"
        component["options"] = lists['atomic_lines']['name']
        component["default"] = component["options"][defaults_lists['atomic_lines']]
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
        component["default"] = False
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
        component["default"] = False
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
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Individual abundance"
        component["options"] = []
        #component["options"].append("1 - H (Hydrogen)") # SPECTRUM has calculations hard coded for H and He
        #component["options"].append("2 - He (Helium)")  # they should not be changed
        component["options"].append("3 - Li (Lithium)")
        component["options"].append("4 - Be (Beryllium)")
        component["options"].append("5 - B (Boron)")
        component["options"].append("6 - C (Carbon)")
        component["options"].append("7 - N (Nitrogen)")
        component["options"].append("8 - O (Oxygen)")
        component["options"].append("9 - F (Fluorine)")
        component["options"].append("10 - Ne (Neon)")
        component["options"].append("11 - Na (Sodium)")
        component["options"].append("12 - Mg (Magnesium)")
        component["options"].append("13 - Al (Aluminium)")
        component["options"].append("14 - Si (Silicon)")
        component["options"].append("15 - P (Phosphorus)")
        component["options"].append("16 - S (Sulfur)")
        component["options"].append("17 - Cl (Chlorine)")
        component["options"].append("18 - Ar (Argon)")
        component["options"].append("19 - K (Potassium)")
        component["options"].append("20 - Ca (Calcium)")
        component["options"].append("21 - Sc (Scandium)")
        component["options"].append("22 - Ti (Titanium)")
        component["options"].append("23 - V (Vanadium)")
        component["options"].append("24 - Cr (Chromium)")
        component["options"].append("25 - Mn (Manganese)")
        component["options"].append("26 - Fe (Iron)") # Determined by MH
        component["options"].append("27 - Co (Cobalt)")
        component["options"].append("28 - Ni (Nickel)")
        component["options"].append("29 - Cu (Copper)")
        component["options"].append("30 - Zn (Zinc)")
        component["options"].append("31 - Ga (Gallium)")
        component["options"].append("32 - Ge (Germanium)")
        component["options"].append("33 - As (Arsenic)")
        component["options"].append("34 - Se (Selenium)")
        component["options"].append("35 - Br (Bromine)")
        component["options"].append("36 - Kr (Krypton)")
        component["options"].append("37 - Rb (Rubidium)")
        component["options"].append("38 - Sr (Strontium)")
        component["options"].append("39 - Y (Yttrium)")
        component["options"].append("40 - Zr (Zirconium)")
        component["options"].append("41 - Nb (Niobium)")
        component["options"].append("42 - Mo (Molybdenum)")
        component["options"].append("43 - Tc (Technetium)")
        component["options"].append("44 - Ru (Ruthenium)")
        component["options"].append("45 - Rh (Rhodium)")
        component["options"].append("46 - Pd (Palladium)")
        component["options"].append("47 - Ag (Silver)")
        component["options"].append("48 - Cd (Cadmium)")
        component["options"].append("49 - In (Indium)")
        component["options"].append("50 - Sn (Tin)")
        component["options"].append("51 - Sb (Antimony)")
        component["options"].append("52 - Te (Tellurium)")
        component["options"].append("53 - I (Iodine)")
        component["options"].append("54 - Xe (Xenon)")
        component["options"].append("55 - Cs (Caesium)")
        component["options"].append("56 - Ba (Barium)")
        component["options"].append("57 - La (Lanthanum)")
        component["options"].append("58 - Ce (Cerium)")
        component["options"].append("59 - Pr (Praseodymium)")
        component["options"].append("60 - Nd (Neodymium)")
        component["options"].append("61 - Pm (Promethium)")
        component["options"].append("62 - Sm (Samarium)")
        component["options"].append("63 - Eu (Europium)")
        component["options"].append("64 - Gd (Gadolinium)")
        component["options"].append("65 - Tb (Terbium)")
        component["options"].append("66 - Dy (Dysprosium)")
        component["options"].append("67 - Ho (Holmium)")
        component["options"].append("68 - Er (Erbium)")
        component["options"].append("69 - Tm (Thulium)")
        component["options"].append("70 - Yb (Ytterbium)")
        component["options"].append("71 - Lu (Lutetium)")
        component["options"].append("72 - Hf (Hafnium)")
        component["options"].append("73 - Ta (Tantalum)")
        component["options"].append("74 - W (Tungsten)")
        component["options"].append("75 - Re (Rhenium)")
        component["options"].append("76 - Os (Osmium)")
        component["options"].append("77 - Ir (Iridium)")
        component["options"].append("78 - Pt (Platinum)")
        component["options"].append("79 - Au (Gold)")
        component["options"].append("80 - Hg (Mercury)")
        component["options"].append("81 - Tl (Thallium)")
        component["options"].append("82 - Pb (Lead)")
        component["options"].append("83 - Bi (Bismuth)")
        component["options"].append("84 - Po (Polonium)")
        component["options"].append("85 - At (Astatine)")
        component["options"].append("86 - Rn (Radon)")
        component["options"].append("87 - Fr (Francium)")
        component["options"].append("88 - Ra (Radium)")
        component["options"].append("89 - Ac (Actinium)")
        component["options"].append("90 - Th (Thorium)")
        component["options"].append("91 - Pa (Protactinium)")
        component["options"].append("92 - U (Uranium)")
        component["options"].append("101 - Md (Mendelevium)")
        component["options"].append("106 - Sg (Seaborgium)")
        component["options"].append("107 - Bh (Bohrium)")
        component["options"].append("108 - Hs (Hassium)")
        component["options"].append("112 - Cn (Copernicium)")
        component["options"].append("113 - Uut (Ununtrium)")
        component["options"].append("114 - Uuq (Ununquadium)")
        component["default"] = component["options"][0]
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "Free individual abundance"
        component["default"] = False
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Maximum number of iterations"
        component["text-type"] = "int" # float, int or str
        component["default"] = "20"
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)


    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


