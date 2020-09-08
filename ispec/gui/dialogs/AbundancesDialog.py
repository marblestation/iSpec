#
#    This file is part of iSpec.
#    Copyright Sergi Blanco-Cuaresma - http://www.blancocuaresma.com/s/
#
#    iSpec is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    iSpec is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with iSpec. If not, see <http://www.gnu.org/licenses/>.
#
import numpy as np
import matplotlib.pyplot as plt
from .CustomDialog import *

class AbundancesDialog(CustomDialog):
    def plot(self, axes, component):
        for i, ax in enumerate(axes):
            ## Draw
            elements = np.unique(self.__linemasks['element'])
            for element in elements:
                flines = self.__linemasks['element'] == element
                element_linemasks = self.__linemasks[flines]
                element_abundances = self.__x_over_h[flines]
                if i == 0:
                    ax.plot(element_linemasks['lower_state_eV'], element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)
                else:
                    ax.plot(element_linemasks['ewr'], element_abundances, linestyle='', marker='o', markersize=5, zorder=1, label=element)

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
        for i in range(len(self.__stats)):
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


    def __init__(self, parent, title, teff, logg, feh, alpha, vmic, lists, default_lists):
        self.__parent = parent
        self.__title = title
        self.__plot = None
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
        component["minvalue"] = 400
        component["maxvalue"] = 55000
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Surface gravity (log g)"
        component["text-type"] = "float" # float, int or str
        component["default"] = logg
        component["minvalue"] = -0.5
        component["maxvalue"] = 5.5
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
        component["type"] = "Entry"
        component["text"] = "Alpha enhancement [alpha/Fe]"
        component["text-type"] = "float" # float, int or str
        component["default"] = alpha
        component["minvalue"] = -2
        component["maxvalue"] = 2
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


