import sys
import Tkinter
import tkMessageBox
from tkSimpleDialog import Dialog

try:
    import ttk
except:
    pass

import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigCanvas, NavigationToolbar2TkAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter

class CustomDialog(Dialog):

    def __init__(self, parent, title, components):

        if not parent:
            parent = Tkinter._default_root

        self.results = None
        self.components = components
        Dialog.__init__(self, parent, title)

    def destroy(self):
        self.components = None
        Dialog.destroy(self)

    def body(self, master):
        first_entry = None
        row = 0
        grid_frame = Tkinter.Frame(master) # Form
        plot_frame = Tkinter.Frame(master) # Plot
        stats_frame = Tkinter.Frame(master) # Listbox
        for i, component in enumerate(self.components):
            row += i
            if component["type"].lower() == "label":
                component["object"] = Tkinter.Label(grid_frame, text=component["text"])
                component["object"].grid(row=row, columnspan=2)
            if component["type"].lower() == "entry":
                Tkinter.Label(grid_frame, text=component["text"]).grid(row=row, sticky=Tkinter.W)
                component["object"] = Tkinter.Entry(grid_frame, justify=Tkinter.RIGHT)
                component["object"].delete(0, Tkinter.END)
                component["object"].insert(0, str(component["default"]))
                component["object"].grid(row=row, column=1)
                if first_entry is None:
                    first_entry = component["object"]
            elif component["type"].lower() == "checkbutton":
                component["variable"] = Tkinter.IntVar() # set to 1 if the button is selected, and 0 otherwise
                component["variable"].set(component["default"])
                component["object"] = Tkinter.Checkbutton(grid_frame, text=component["text"], variable=component["variable"])
                component["object"].grid(row=row, columnspan=2, sticky=Tkinter.E)
            #elif component["type"].lower() == "combobox":
            elif component["type"].lower() == "optionmenu":
                if len(component["options"]) <= 0:
                    raise Exception("Options needed!")
                Tkinter.Label(grid_frame, text=component["text"]).grid(row=row, sticky=Tkinter.W)
                component["variable"] = Tkinter.StringVar()

                if "ttk" in sys.modules.keys():
                    # It supports better long lists of options
                    component["variable"].set(component["default"])
                    component["object"] = ttk.Combobox(grid_frame, textvariable=component["variable"], state='readonly')
                    component["object"]['values'] = tuple(component["options"])
                    component["object"].grid(row=row, column=1, sticky=Tkinter.W)
                else:
                    Tkinter.Label(grid_frame, text=component["text"]).grid(row=row, sticky=Tkinter.W)
                    component["variable"] = Tkinter.StringVar()
                    component["variable"].set(component["default"])
                    component["object"] = apply(Tkinter.OptionMenu, (grid_frame, component["variable"]) + tuple(component["options"]))
                    component["object"].grid(row=row, column=1, sticky=Tkinter.W)
            elif component["type"].lower() == "radiobutton":
                Tkinter.Label(grid_frame, text=component["text"]).grid(row=row, sticky=Tkinter.W)
                component["variable"] = Tkinter.StringVar()
                objects = []
                for j, text in enumerate(component["options"]):
                    objects.append(Tkinter.Radiobutton(grid_frame, text=text, variable=component["variable"], value=text))
                    objects[-1].grid(row=row, column=1, sticky=Tkinter.W)
                    row += 1
                component["variable"].set(component["default"])
            elif component["type"].lower() == "listbox" and component["options"] is not None and len(component["options"]) > 0:
                stats_scrollbar = Tkinter.Scrollbar(stats_frame, orient=Tkinter.VERTICAL)
                component["object"] = Tkinter.Listbox(stats_frame, height=5, yscrollcommand=stats_scrollbar.set, font=('courier',10,'normal'), selectmode=Tkinter.EXTENDED)
                stats_scrollbar.pack(side=Tkinter.RIGHT, fill=Tkinter.Y)
                component["object"].delete(0, Tkinter.END) # Delete all
                for text in component["options"]:
                    component["object"].insert(Tkinter.END, text)
                component["object"].pack(fill=Tkinter.BOTH, expand=1)
            elif component["type"].lower() == "plot" and component["function"] is not None:
                # Create the mpl Figure and FigCanvas objects.
                # 5x4 inches, 100 dots-per-inch
                #
                self.dpi = 100
                self.fig = Figure((6.0, 3.0), dpi=self.dpi)
                self.canvas = FigCanvas(self.fig, master=plot_frame)
                self.canvas.show()
                #self.canvas.get_tk_widget().grid(row = 0, column = 0)  # or .grid(row = 0)

                ### Since we have only one plot, we can use add_axes
                ### instead of add_subplot, but then the subplot
                ### configuration tool in the navigation toolbar wouldn't
                ### work.
                ###
                if component["axes"] == 1:
                    self.axes = self.fig.add_subplot(1, 1, 1)
                    #### Avoid using special notation that are not easy to understand in axis for big zoom
                    myyfmt = ScalarFormatter(useOffset=False)
                    self.axes.get_xaxis().set_major_formatter(myyfmt)
                    self.axes.get_yaxis().set_major_formatter(myyfmt)

                    self.toolbar = NavigationToolbar(self.canvas, plot_frame)
                    self.toolbar.update()

                    self.canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
                elif component["axes"] > 1:
                    self.axes = []
                    for i in xrange(component["axes"]):
                        self.axes.append(self.fig.add_subplot(2, 1, i+1))
                        #### Avoid using special notation that are not easy to understand in axis for big zoom
                        myyfmt = ScalarFormatter(useOffset=False)
                        self.axes[i].get_xaxis().set_major_formatter(myyfmt)
                        self.axes[i].get_yaxis().set_major_formatter(myyfmt)

                    self.toolbar = NavigationToolbar(self.canvas, plot_frame)
                    self.toolbar.update()

                    self.canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

                # Plotting function
                apply(component["function"], (self.axes, component))

        plot_frame.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        stats_frame.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        grid_frame.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

        # Which should receive focus:
        return first_entry

    def validate(self):

        for component in self.components:
            if component["type"] == "Entry":
                try:
                    #value = float(component["object"].get())
                    if component["text-type"] == "str":
                        value = eval(component["text-type"] + "('%s')" % component["object"].get())
                    else:
                        value = eval(component["text-type"] + "(%s)" % component["object"].get())
                except Exception:
                    tkMessageBox.showwarning(
                        "Illegal value",
                        "Illegal value for " + component["text"] + "\nPlease try again",
                        parent = self
                    )
                    return 0
                if component["minvalue"] is not None and value < component["minvalue"]:
                    tkMessageBox.showwarning(
                        "Too small",
                        "'%s': "
                        "The allowed minimum value is %s. "
                        "Please try again." % (component["text"], component["minvalue"]),
                        parent = self
                    )
                    return 0
                if component["maxvalue"] is not None and value > component["maxvalue"]:
                    tkMessageBox.showwarning(
                        "Too big",
                        "'%s': "
                        "The allowed maximum value is %s. "
                        "Please try again." % (component["text"], component["maxvalue"]),
                        parent = self
                    )
                    return 0
        return 1

    def apply(self):
        self.results = {}
        for component in self.components:
            # Save results and modify defaults
            if component["type"] == "Entry":
                if component["text-type"] == "str":
                    value = eval(component["text-type"] + "('%s')" % component["object"].get())
                else:
                    value = eval(component["text-type"] + "(%s)" % component["object"].get())
                self.results[component["text"]] = value
                component["default"] = str(value)
            elif component["type"] in ["Checkbutton", "OptionMenu", "Radiobutton"]:
                self.results[component["text"]] = component["variable"].get()
                component["default"] = component["variable"].get()


