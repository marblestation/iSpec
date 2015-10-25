import numpy as np
import matplotlib.pyplot as plt
from CustomDialog import *

class AbundancesDialog(CustomDialog):
    def plot(self, axes, component):
        ## Draw
        elements = np.unique(self.__linemasks['element'])
        for element in elements:
            flines = self.__linemasks['element'] == element
            element_linemasks = self.__linemasks[flines]
            element_abundances = self.__x_over_h[flines]
            #axes.plot(element_linemasks['VALD_wave_peak'], element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)
            axes.plot(element_linemasks['lower_state_eV'], element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)

        leg = axes.legend(loc='upper right', shadow=False)
        ltext  = leg.get_texts()
        plt.setp(ltext, fontsize='8')

        axes.grid(True, which="both")
        axes.set_title("Abundances", fontsize="10")
        #axes.set_xlabel("wavelength (nm)", fontsize="10")
        axes.set_xlabel("lower state (eV)", fontsize="10")
        axes.set_ylabel("abundance (dex)", fontsize="10")
        fig = axes.get_figure()
        fig.set_tight_layout(True)

    def register(self, linemasks, x_over_h, x_over_fe):
        self.__x_over_h = x_over_h
        self.__x_over_fe = x_over_fe
        self.__linemasks = linemasks

        # We have data, we can assign the plotting function
        self.__components[0]["function"] = self.plot

        teff = float(self.__components[5]["default"])
        logg = float(self.__components[6]["default"])
        feh = float(self.__components[7]["default"])
        vmic = float(self.__components[8]["default"])

        ## Stats
        for i in xrange(len(self.__stats)):
            self.__stats.pop()
        self.__stats.append("%-50s: %10.2f" % ("Effective temperature (k)", np.round(teff, 1)))
        self.__stats.append("%-50s: %10.2f" % ("Surface gravity (log g)", np.round(logg, 2)))
        self.__stats.append("%-50s: %10.2f" % ("Metallicity [Fe/H]", np.round(feh, 2)))
        self.__stats.append("%-50s: %10.2f" % ("Microturbulence velocity (km/s)", np.round(vmic, 2)))
        self.__stats.append("%-50s: %10.2f" % ("Total number of lines", np.round(len(x_over_h), 2)))

        elements = np.unique(linemasks['element'])
        for element in elements:
            flines = linemasks['element'] == element
            element_linemasks = linemasks[flines]
            element_abundances_over_h = x_over_h[flines]
            element_abundances_over_fe = x_over_fe[flines]
            self.__stats.append("%-50s: %10.2f" % ( element + " median abundance in [X/H] (dex)", np.round(np.nanmedian(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " mean abundance in [X/H] (dex)", np.round(np.nanmean(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " standard deviation in [X/H] (dex)", np.round(np.nanstd(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " median abundance in [X/Fe] (dex)", np.round(np.nanmedian(element_abundances_over_fe), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " mean abundance in [X/Fe] (dex)", np.round(np.nanmean(element_abundances_over_fe), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " standard deviation in [X/Fe] (dex)", np.round(np.nanstd(element_abundances_over_fe), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " lines number", np.round(len(element_abundances_over_h), 0)))


    def __init__(self, parent, title, teff, logg, feh, vmic, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__plot = None
        self.__stats = []
        self.__components = []
        component = {}
        component["type"] = "Plot"
        component["function"] = self.__plot
        self.__components.append(component)
        component = {}
        component["type"] = "Listbox"
        component["options"] = self.__stats
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Code"
        component["options"] = lists['ew_code']
        component["default"] = component["options"][default_lists['ew_code']]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Model atmosphere"
        component["options"] = lists['atmospheres']['name']
        component["default"] = component["options"][default_lists['atmospheres']]
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Solar abundances"
        component["options"] = lists['abundances']['name']
        component["default"] = component["options"][default_lists['abundances']]
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

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


