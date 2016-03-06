import numpy as np
from CustomDialog import *

class VelocityProfileDialog(CustomDialog):
    def plot(self, axes, component):
        ## Draw
        axes.plot(self.__xcoord, self.__fluxes, lw=1, color='b', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)
        axes.fill_between(self.__xcoord, self.__fluxes+self.__errors, self.__fluxes-self.__errors, color='#CCCCCC')

        if self.__velocity_step >= 0.1:
            xcoord_mod = np.arange(np.min(self.__xcoord), np.max(self.__xcoord), 0.1)
            for model in self.__models:
                axes.plot(xcoord_mod, model(xcoord_mod), lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)
        else:
            for model in self.__models:
                axes.plot(self.__xcoord, model(self.__xcoord), lw=1, color='r', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='r', zorder=2)

        axes.grid(True, which="both")
        axes.set_title("Profile", fontsize="10")
        axes.set_xlabel("velocity (km/s)", fontsize="10")
        axes.set_ylabel("relative intensity", fontsize="10")
        fig = axes.get_figure()
        fig.set_tight_layout(True)

    def register(self, ccf, models, telluric_fwhm=None):
        self.__xcoord = ccf['x']
        self.__fluxes = ccf['y']
        self.__errors = ccf['err']
        self.__models = models
        # The default values have been previously updated with the user selected values
        self.__velocity_step = float(self.__components[4]["default"])
        if len(self.__components) == 9 and self.__components[8]["text"] == "Cross-correlate with":
            self.__template = self.__components[8]["default"]
        else:
            self.__template = None

        # We have data, we can assign the plotting function
        self.__components[0]["function"] = self.plot

        ## Stats
        for i in xrange(len(self.__stats)):
            self.__stats.pop()
        if self.__template is not None:
            self.__stats.append("%-50s: %s" % ("Template", str(self.__template)))
        for i, model in enumerate(models):
            self.__stats.append("%-50s: %10.2f" % ("Mean (km/s)", np.round(model.mu(), 2)))
            self.__stats.append("%-50s: %10.2f" % ("Error (+/- km/s)", np.round(model.emu(), 2)))
            self.__stats.append("%-50s: %10.2f" % ("Baseline", np.round(model.baseline(), 2)))
            self.__stats.append("%-50s: %10.2f" % ("A (rel. intensity)", np.round(model.A(), 2)))
            self.__stats.append("%-50s: %10.2f" % ("Sigma (km/s)", np.round(model.sig(), 2)))

            try:
                # If model is VoigtModel
                self.__stats.append("%-50s: %10.2f" % ("Gamma", np.round(model.gamma(), 2)))
            except AttributeError:
                # model is GaussianModel
                pass

            fwhm = model.fwhm()[0] # km/s (because xcoord is already velocity)
            self.__stats.append("%-50s: %10.2f" % ("Measured FWHM (km/s)", np.round(fwhm, 2)))
            if telluric_fwhm != 0:
                self.__stats.append("%-50s: %10.2f" % ("Theoretical telluric FWHM (km/s)", np.round(telluric_fwhm, 2)))
                self.__stats.append("%-50s: %10.2f" % ("Corrected FWHM (km/s)", np.round(fwhm - telluric_fwhm, 2)))
            c = 299792458.0 # m/s
            if (fwhm - telluric_fwhm > 0.0):
                R = np.int(c/(1000.0*(fwhm - telluric_fwhm)))
                self.__stats.append("%-50s: %10.i" % ("Estimated resolving power (R)", np.round(R)))

            self.__stats.append("%-50s: %10.5f" % ("RMS", np.round(model.rms, 5)))
            self.__stats.append("%-50s  %s" % ("------------------------------", "----------"))

    def __init__(self, parent, title, rv_upper_limit, rv_lower_limit, rv_step, templates, masks, mask_size=2.0, mask_depth=0.01):
        self.__parent = parent
        self.__title = title
        self.__plot = None
        self.__stats = []
        self.__components = []
        component = {}
        component["type"] = "Plot"
        component["function"] = self.__plot
        component["axes"] = 1
        self.__components.append(component)
        component = {}
        component["type"] = "Listbox"
        component["options"] = self.__stats
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity lower limit (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = rv_lower_limit
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity upper limit (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = rv_upper_limit
        component["minvalue"] = -np.inf
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Velocity steps (km/s)"
        component["text-type"] = "float" # float, int or str
        component["default"] = rv_step
        component["minvalue"] = 0.001
        component["maxvalue"] = np.inf
        self.__components.append(component)
        component = {}
        component["type"] = "OptionMenu"
        component["text"] = "Fitting model"
        component["options"] = ['2nd order polynomial + auto fit', '2nd order polynomial + gaussian fit', '2nd order polynomial + voigt fit']
        component["default"] = component["options"][1]
        self.__components.append(component)
        component = {}
        component["type"] = "Entry"
        component["text"] = "Peak probability"
        component["text-type"] = "float" # float, int or str
        component["default"] = 0.75
        component["minvalue"] = 0.0
        component["maxvalue"] = 1.0
        self.__components.append(component)
        component = {}
        component["type"] = "Checkbutton"
        component["text"] = "CCF in Fourier space"
        component["default"] = False
        self.__components.append(component)
        if len(templates) > 0:
            component = {}
            component["type"] = "OptionMenu"
            component["text"] = "Cross-correlate with"
            component["options"] = templates
            component["default"] = templates[0]
            self.__components.append(component)
        else:
            if len(masks) > 1:
                component = {}
                component["type"] = "OptionMenu"
                component["text"] = "Mask linelist"
                component["options"] = masks
                component["default"] = masks[0]
                self.__components.append(component)
            component = {}
            component["type"] = "Entry"
            component["text"] = "Mask size (km/s)"
            component["text-type"] = "float" # float, int or str
            component["default"] = mask_size
            component["minvalue"] = 0.001
            component["maxvalue"] = np.inf
            self.__components.append(component)
            component = {}
            component["type"] = "Entry"
            component["text"] = "Minimum depth"
            component["text-type"] = "float" # float, int or str
            component["default"] = mask_depth
            component["minvalue"] = 0.001
            component["maxvalue"] = 1.0
            self.__components.append(component)

    def show(self, updated_templates=None):
        self.results = None
        if updated_templates is not None:
            self.__components[8]["options"] = updated_templates
            # Validate that the default value (previous user selected value) exists in the new template list
            default_ok = False
            for template in updated_templates:
                if template == self.__components[8]["default"]:
                    default_ok = True
            if not default_ok:
                self.__components[8]["default"] = updated_templates[0]
        CustomDialog.__init__(self, self.__parent, self.__title, self.__components)

