import numpy as np
import matplotlib.pyplot as plt
from CustomDialog import *

class SolverEWDialog(CustomDialog):
    def plot(self, axes, component):
        for i, ax in enumerate(axes):
            ## Draw
            for j, element in enumerate(["Fe 1", "Fe 2"]):
                element_linemasks = self.__linemasks[self.__linemasks['element'] == element]
                element_abundances = self.__x_over_h[self.__selected_x_over_h[j]]
                if len(element_abundances) == 0:
                    continue
                if i == 0:
                    m = self.__fitted_lines_params[0]
                    c = self.__fitted_lines_params[1]
                    x = element_linemasks['lower_state_eV']
                else:
                    m = self.__fitted_lines_params[2]
                    c = self.__fitted_lines_params[3]
                    x = element_linemasks['ewr']
                ax.plot(x, element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)
                ax.plot(x, m*x + c, color="red")

            leg = ax.legend(loc='upper right', shadow=False, numpoints=1)
            ltext  = leg.get_texts()
            plt.setp(ltext, fontsize='8')

            ax.grid(True, which="both")
            if i == 0:
                ax.set_xlabel("lower state (eV)", fontsize="10")
            else:
                ax.set_xlabel("reduced equivalent width", fontsize="10")
            ax.set_ylabel("[X/H]", fontsize="10")
            ax.tick_params(axis='y', labelsize=8)
        fig = ax.get_figure()
        #fig.set_tight_layout(True)
        fig.subplots_adjust(hspace = 0.5, bottom=0.2)

    def register(self, linemasks, params, x_over_h, selected_x_over_h, fitted_lines_params):
        self.__linemasks = linemasks
        self.__params = params
        self.__x_over_h = x_over_h
        self.__selected_x_over_h = selected_x_over_h
        self.__fitted_lines_params = fitted_lines_params

        # We have data, we can assign the plotting function
        self.__components[0]["function"] = self.plot

        teff = params['teff']
        logg = params['logg']
        feh = params['MH']
        vmic = params['vmic']

        ## Stats
        for i in xrange(len(self.__stats)):
            self.__stats.pop()
        self.__stats.append("%-50s: %10.2f" % ("Effective temperature (k)", np.round(teff, 1)))
        self.__stats.append("%-50s: %10.2f" % ("Surface gravity (log g)", np.round(logg, 2)))
        self.__stats.append("%-50s: %10.2f" % ("Metallicity [Fe/H]", np.round(feh, 2)))
        self.__stats.append("%-50s: %10.2f" % ("Microturbulence velocity (km/s)", np.round(vmic, 2)))

        self.__stats.append("%-50s: %10.4f" % ("Excitation potential slope", np.round(self.__fitted_lines_params[0], 2)))
        self.__stats.append("%-50s: %10.4f" % ("Reduced equivalent width slope", np.round(self.__fitted_lines_params[2], 2)))
        diff = np.nanmedian(self.__x_over_h[self.__selected_x_over_h[0]]) - np.nanmedian(self.__x_over_h[self.__selected_x_over_h[1]])
        self.__stats.append("%-50s: %10.2f" % ("Fe I - Fe II abundance difference", np.round(diff, 2)))

        self.__stats.append("%-50s: %10.0f" % ("Total number of lines", np.round(len(x_over_h), 2)))

        for i, element in enumerate(["Fe 1", "Fe 2"]):
            element_abundances_over_h = x_over_h[selected_x_over_h[i]]
            if len(element_abundances_over_h) == 0:
                continue
            self.__stats.append("%-50s: %10.2f" % ( element + " median abundance in [X/H] (dex)", np.round(np.nanmedian(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " mean abundance in [X/H] (dex)", np.round(np.nanmean(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.2f" % ( element + " standard deviation in [X/H] (dex)", np.round(np.nanstd(element_abundances_over_h), 2)))
            self.__stats.append("%-50s: %10.0f" % ( element + " lines number", np.round(len(element_abundances_over_h), 0)))


    def __init__(self, parent, title, teff, logg, feh, vmic, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__plot = None
        self.__params = None
        self.__x_over_h = None
        self.__selected_x_over_h = None
        self.__fitted_lines_params = None
        self.__stats = []
        self.__components = []
        component = {}
        component["type"] = "Plot"
        component["function"] = self.__plot
        component["axes"] = 2
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
        component["text"] = "Maximum number of iterations"
        component["text-type"] = "int" # float, int or str
        component["default"] = "10"
        component["minvalue"] = 0
        component["maxvalue"] = np.inf
        self.__components.append(component)

    def show(self):
        self.results = None
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)


