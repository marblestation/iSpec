#!/usr/bin/env python
import Tkinter
import tkMessageBox
import tkFileDialog
import tkSimpleDialog
#import matplotlib
#matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigCanvas, NavigationToolbar2TkAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt

import os
import sys

import random
import numpy as np
import threading
import logging


import sve
from dialogs import *
from CustomizableRegion import *
from Spectrum import *
from Meter import *
from StatusBar import *


try:
    from SAMPManager import *
except: pass

def resource_path(relative):
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
        # Since we are inside "sve/gui/", we go up two levels to find "input/"
        basedir = os.path.dirname(basedir[:-1])
        basedir = os.path.dirname(basedir[:-1])
    return os.path.join(basedir, relative)


class SVEBaseApp(Tkinter.Tk):

    def __init_attributes__(self):
        self.velocity_telluric_lower_limit = -100 # km/s
        self.velocity_telluric_upper_limit = 100 # km/s
        self.velocity_telluric_step = 0.5 # km/s
        self.linelist_telluric = None
        self.velocity_atomic_lower_limit = -200 # km/s
        self.velocity_atomic_upper_limit = 200 # km/s
        self.velocity_atomic_step = 1.0 # km/s
        self.linelist_atomic = None
        self.velocity_template_lower_limit = -200 # km/s
        self.velocity_template_upper_limit = 200 # km/s
        self.velocity_template_step = 1.0 # km/s
        self.modeled_layers_pack = {} # Synthesize spectrum (atmospheric models)
        self.find_continuum_regions_wave_step = 0.05
        self.find_continuum_regions_sigma = 0.001
        self.find_continuum_regions_max_continuum_diff = 1.0
        self.find_lines_min_depth = 0.05 # (% of the continuum)
        self.find_lines_max_depth = 1.00 # (% of the continuum)
        self.show_errors = False

        self.linelist_SPECTRUM = {}
        self.abundances_SPECTRUM = {}

        # Barycentric velocity determination (from date, time, right ascension and declination)
        self.barycentric_vel = 0.0 # km/s

        # Operations
        self.operation_waveobs = "waveobs"
        self.operation_flux = "flux"
        self.operation_err = "err"
        self.safe_operations_description = []
        self.safe_operations = {}
        self.safe_operations_description.append("sin(x)  ::  Trigonometric sine, element-wise.")
        self.safe_operations['sin'] = np.sin
        self.safe_operations_description.append("cos(x)  ::  Cosine elementwise.")
        self.safe_operations['cos'] = np.cos
        self.safe_operations_description.append("tan(x)  ::  Compute tangent element-wise.")
        self.safe_operations['tan'] = np.tan
        self.safe_operations_description.append("arcsin(x)  ::  Inverse sine, element-wise.")
        self.safe_operations['arcsin'] = np.arcsin
        self.safe_operations_description.append("arccos(x)  ::  Trigonometric inverse cosine, element-wise.")
        self.safe_operations['arccos'] = np.arccos
        self.safe_operations_description.append("arctan(x)  ::  Trigonometric inverse tangent, element-wise.")
        self.safe_operations['arctan'] = np.arctan
        self.safe_operations_description.append("arctan2(x1, x2)  ::  Arc tangent of x1/x2 choosing the correct quadrant")
        self.safe_operations['arctan2'] = np.arctan2

        self.safe_operations_description.append("sinh(x)  ::  Hyperbolic sine, element-wise.")
        self.safe_operations['sinh'] = np.sinh
        self.safe_operations_description.append("cosh(x)  ::  Hyperbolic cosine, element-wise.")
        self.safe_operations['cosh'] = np.cosh
        self.safe_operations_description.append("tanh(x)  ::  Compute hyperbolic tangent element-wise.")
        self.safe_operations['tanh'] = np.tanh
        self.safe_operations_description.append("arcsinh(x)  ::  Inverse hyperbolic sine elementwise.")
        self.safe_operations['arcsinh'] = np.arcsinh
        self.safe_operations_description.append("arccosh(x)  ::  Inverse hyperbolic cosine, elementwise.")
        self.safe_operations['arccosh'] = np.arccosh
        self.safe_operations_description.append("arctanh(x)  ::  Inverse hyperbolic tangent elementwise.")
        self.safe_operations['arctanh'] = np.arctanh

        self.safe_operations_description.append("around(a[, dec])  ::  Evenly round to the given number of decimals.")
        self.safe_operations['around'] = np.around
        self.safe_operations_description.append("floor(x)  ::  Return the floor of the input, element-wise.")
        self.safe_operations['floor'] = np.floor
        self.safe_operations_description.append("ceil(x)  ::  Return the ceiling of the input, element-wise.")
        self.safe_operations['ceil'] = np.ceil

        self.safe_operations_description.append("exp(x)  ::  Calculate the exponential of all elements in the input array.")
        self.safe_operations['exp'] = np.exp
        self.safe_operations_description.append("log(x)  ::  Natural logarithm, element-wise.")
        self.safe_operations['log'] = np.log
        self.safe_operations_description.append("log10(x)  ::  Return the base 10 logarithm of the input array, element-wise.")
        self.safe_operations['log10'] = np.log10
        self.safe_operations_description.append("log2(x)  ::  Base-2 logarithm of x.")
        self.safe_operations['log2'] = np.log2

        self.safe_operations_description.append("sqrt(x)  ::  Return the positive square-root of an array, element-wise.")
        self.safe_operations['sqrt'] = np.sqrt
        self.safe_operations_description.append("absolute(x)  ::  Compute the absolute values elementwise.")
        self.safe_operations['absolute'] = np.absolute

        self.safe_operations_description.append("add(x1, x2)  ::  Add arguments element-wise.")
        self.safe_operations['add'] = np.add
        self.safe_operations_description.append("multiply(x1, x2)  ::  Multiply arguments element-wise.")
        self.safe_operations['multiply'] = np.multiply
        self.safe_operations_description.append("divide(x1, x2)  ::  Divide arguments element-wise.")
        self.safe_operations['divide'] = np.divide
        self.safe_operations_description.append("power(x1, x2)  ::  Array elements raised to 'x2' powers.")
        self.safe_operations['power'] = np.power
        self.safe_operations_description.append("subtract(x1, x2)  ::  Subtract arguments, element-wise.")
        self.safe_operations['subtract'] = np.subtract
        self.safe_operations_description.append("mod(x1, x2)  ::  Return element-wise remainder of division.")
        self.safe_operations['mod'] = np.mod

        self.safe_operations_description.append("pi  ::  Pi value.")
        self.safe_operations['pi'] = np.pi
        self.safe_operations_description.append("e  ::  Euler's number.")
        self.safe_operations['e'] = np.e

    def __init__(self, spectra, regions, filenames):
        Tkinter.Tk.__init__(self)
        self.protocol('WM_DELETE_WINDOW', self.on_close)

        # Window icon
        img = Tkinter.PhotoImage(file=resource_path("images/SVE.gif"))
        self.tk.call('wm', 'iconphoto', self._w, img)
        #self.iconbitmap(bitmap="@"+resource_path("images/SVE.xbm")) # Black and white


        self.__init_attributes__()

        self.spectra_colors = ('#0000FF', '#A52A2A', '#A020F0', '#34764A', '#000000', '#90EE90', '#FFA500', '#1E90FF',   '#FFC0CB', '#7F7F7F', '#00FF00',)

        self.dupiclated_name_separator = "#"
        self.region_widgets = {} # It is important not to lose the reference to the regions or the garbage
        self.region_widgets['continuum'] = []
        self.region_widgets['lines'] = []
        self.region_widgets['segments'] = []

        self.action = "Stats"
        self.elements = "continuum"
        self.subelements = None
        self.lock = threading.Lock() # Used by CustomRegion
        self.not_saved = {}
        self.not_saved["continuum"] = False
        self.not_saved["lines"] = False
        self.not_saved["segments"] = False
        self.show_errors = False

        self.operation_in_progress = False
        self.programmed_flash_status = None
        self.spectrum_function_items = []
        self.dialog = {}

        self.spectra = []
        self.active_spectrum = None
        self.active_spectrum_history = [] # History to be able to come back to the previous active when the current is closed
        for i in np.arange(len(spectra)):
            path = filenames["spectra"][i]
            name = path.split('/')[-1]
            name = self.get_name(name) # If it already exists, add a suffix
            color = self.get_color()
            self.active_spectrum = Spectrum(spectra[i], name, path=path, color=color)
            self.spectra.append(self.active_spectrum)
            self.active_spectrum_history.append(self.active_spectrum)


        if regions == None:
            continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
            lines = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
            segments = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])

            self.regions = {}
            self.regions['continuum'] = continuum
            self.regions['lines'] = lines
            self.regions['segments'] = segments
        else:
            self.regions = regions

        if filenames == None:
            self.filenames = {}
            self.filenames['continuum'] = None
            self.filenames['lines'] = None
            self.filenames['segments'] = None
        else:
            self.filenames = {}
            self.filenames['continuum'] = filenames['continuum']
            self.filenames['lines'] = filenames['lines']
            self.filenames['segments'] = filenames['segments']


        try:
            self.samp_manager = SAMPManager(self.on_receive_spectrum, check_connection_period=2)
        except:
            self.samp_manager = None

        self.create_main_window()
        self.draw_figure_for_first_time()

        # In MacOSX the window does not get the focus after execution from terminal
        #if sys.platform == "darwin":
            ## - Iconifying and deiconifying makes the user see some changes in the screen at least
            #self.iconify()
            #self.update()
            #self.deiconify()
            #self.tkraise()
            #self.focus_set()
            #self.focus_force()


    def create_main_window(self):
        self.create_window()
        self.create_menu()
        self.create_plotting_area()
        self.create_controls()
        self.create_statusbar()

    def create_window(self ):
        self.frame = Tkinter.Frame(self)
        self.frame.pack(fill=Tkinter.BOTH, expand=1)
        self.wm_title("SVE")

    def create_menu(self):
        # create a menu
        menu = Tkinter.Menu(self)
        self.config(menu=menu)

        filemenu = Tkinter.Menu(menu)
        menu.add_cascade(label="Files", menu=filemenu)
        filemenu.add_command(label="Open spectra", command=self.on_open_spectra)
        filemenu.add_command(label="Open continuum regions", command=self.on_open_continuum)
        filemenu.add_command(label="Open line regions", command=self.on_open_lines)
        filemenu.add_command(label="Open segment regions", command=self.on_open_segments)
        filemenu.add_separator()
        filemenu.add_command(label="Save plot image as...", command=self.on_save_plot)
        filemenu.add_command(label="Save spectrum as...", command=self.on_save_spectrum)
        self.spectrum_function_items.append((filemenu, filemenu.entrycget(Tkinter.END, "label")))
        filemenu.add_command(label="Save continuum regions as...", command=self.on_save_continuum_regions)
        filemenu.add_command(label="Save line regions as...", command=self.on_save_line_regions)
        filemenu.add_command(label="Save segment as...", command=self.on_save_segments)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.on_close)

        operationmenu = Tkinter.Menu(menu)
        menu.add_cascade(label="Operations", menu=operationmenu)
        operationmenu.add_command(label="Fit continuum", command=self.on_fit_continuum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Fit lines", command=self.on_fit_lines)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_separator()

        clearmenu = Tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Clear...", menu=clearmenu)
        clearmenu.add_command(label="Fitted continuum", command=self.on_remove_fitted_continuum)
        clearmenu.add_command(label="Fitted lines", command=self.on_remove_fitted_lines)
        clearmenu.add_command(label="Continuum regions", command=self.on_remove_continuum_regions)
        clearmenu.add_command(label="Line regions", command=self.on_remove_line_masks)
        clearmenu.add_command(label="Segments", command=self.on_remove_segments)

        operationmenu.add_separator()
        operationmenu.add_command(label="Find continuum regions", command=self.on_find_continuum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Find line masks", command=self.on_find_lines)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))

        operationmenu.add_separator()
        velocitymenu = Tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Determine velocity relative to...", menu=velocitymenu)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))

        velocitymenu.add_command(label="Atomic line mask (radial velocity)", command=self.on_determine_velocity_atomic)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))
        velocitymenu.add_command(label="Telluric line mask  (barycentric velocity)", command=self.on_determine_velocity_telluric)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))
        velocitymenu.add_command(label="Template", command=self.on_determine_velocity_template)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))


        velocitymenu = Tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Correct velocity relative to...", menu=velocitymenu)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))

        velocitymenu.add_command(label="Atomic line mask (radial velocity)", command=self.on_correct_velocity_atomic)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))
        velocitymenu.add_command(label="Telluric line mask  (barycentric velocity)", command=self.on_correct_velocity_telluric)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))
        velocitymenu.add_command(label="Template", command=self.on_correct_velocity_template)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(Tkinter.END, "label")))


        operationmenu.add_command(label="Calculate barycentric velocity", command=self.on_determine_barycentric_vel)
        operationmenu.add_separator()
        operationmenu.add_command(label="Estimate SNR", command=self.on_estimate_snr)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Estimate errors based on SNR", command=self.on_estimate_errors)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_separator()
        operationmenu.add_command(label="Degrade resolution", command=self.on_degrade_resolution)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Continuum normalization", command=self.on_continuum_normalization)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Clean fluxes and errors", command=self.on_clean_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Clean telluric regions", command=self.on_clean_tellurics)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Wavelength range reduction", command=self.on_cut_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Apply mathematical expression", command=self.on_operate_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Resample spectrum", command=self.on_resample_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))
        operationmenu.add_command(label="Combine all spectra", command=self.on_combine_spectra)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(Tkinter.END, "label")))


        if "determine_abundances" in dir(sve):
            parametersmenu = Tkinter.Menu(menu)
            menu.add_cascade(label="Parameters", menu=parametersmenu)
            #parametersmenu.add_command(label="Find astrophysical parameters", command=self.on_determine_parameters, state=Tkinter.DISABLED)
            #self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(Tkinter.END, "label")))
            parametersmenu.add_command(label="Determine abundances with fitted lines", command=self.on_determine_abundances)
            self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(Tkinter.END, "label")))


        self.menu_active_spectrum_num = Tkinter.IntVar()
        self.menu_active_spectrum_num.set('1')

        self.menu_active_spectrum = Tkinter.Menu(menu)
        menu.add_cascade(label="Spectra", menu=self.menu_active_spectrum)

        self.menu_active_spectrum.add_command(label="Duplicate spectrum", command=self.on_duplicate_spectrum)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(Tkinter.END, "label")))
        self.menu_active_spectrum.add_command(label="Close spectrum", command=self.on_close_spectrum)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(Tkinter.END, "label")))
        self.menu_active_spectrum.add_command(label="Close all spectra", command=self.on_close_all_spectra)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(Tkinter.END, "label")))
        self.menu_active_spectrum.add_separator()

        if "generate_spectrum" in dir(sve):
                self.menu_active_spectrum.add_command(label="Synthesize spectrum", command=self.on_synthesize)

        if self.samp_manager != None:
            self.menu_active_spectrum.add_command(label="Send spectrum to...", command=self.on_send_spectrum)
            self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(Tkinter.END, "label")))
        self.show_errors = Tkinter.BooleanVar()
        self.menu_active_spectrum.add_checkbutton(label="Show errors in plot", onvalue=True, offvalue=False, variable=self.show_errors, command=self.on_show_errors)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(Tkinter.END, "label")))
        self.menu_active_spectrum.add_separator()

        helpmenu = Tkinter.Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)

        helpmenu.add_command(label="License", command=self.on_license)
        helpmenu.add_command(label="About...", command=self.on_about)
        self.update_menu_active_spectrum()

    def create_plotting_area(self):
        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.plot_frame = Tkinter.Frame(self.frame)
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.fig, master=self.plot_frame)
        self.canvas.show()

        # Since we have only one plot, we can use add_axes
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(1, 1, 1)
        # Avoid using special notation that are not easy to understand in axis for big zoom
        myyfmt = ScalarFormatter(useOffset=False)
    	self.axes.get_xaxis().set_major_formatter(myyfmt)
    	self.axes.get_yaxis().set_major_formatter(myyfmt)

        # Make space for the legend
        box = self.axes.get_position()
        self.axes.set_position([box.x0, box.y0, box.width * 0.80, box.height])

        #self.canvas.mpl_connect('button_release_event', self.on_release)
        #self.canvas.mpl_connect('motion_notify_event', self.on_motion)

        self.toolbar = NavigationToolbar(self.canvas, self.plot_frame)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)
        self.plot_frame.pack(side=Tkinter.TOP, fill=Tkinter.BOTH, expand=1)

    def create_controls(self):
        # controls
        self.control_frame = Tkinter.Frame(self.frame)

        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

        label = Tkinter.Label(self.control_frame, text="Action:")
        label.pack(side=Tkinter.LEFT)
        self.action_entry = Tkinter.StringVar()
        self.action_entry.set("Stats")
        self.action_buttons = {}
        for text in ["Stats", "Create", "Modify", "Remove"]:
            self.action_buttons[text] = Tkinter.Radiobutton(self.control_frame, text=text, variable=self.action_entry, value=text, command=self.on_action_change)
            self.action_buttons[text].pack(side=Tkinter.LEFT)

        label = Tkinter.Label(self.control_frame, text=" || Element:")
        label.pack(side=Tkinter.LEFT)
        self.elements_entry = Tkinter.StringVar()
        self.elements_entry.set("Continuum")
        self.elements_buttons = {}
        for text in ["Continuum", "Lines", "Line marks", "Segments"]:
            if text == "Line marks":
                self.elements_buttons[text] = Tkinter.Radiobutton(self.control_frame, text=text, variable=self.elements_entry, value=text, command=self.on_element_change, state=Tkinter.DISABLED)
            else:
                self.elements_buttons[text] = Tkinter.Radiobutton(self.control_frame, text=text, variable=self.elements_entry, value=text, command=self.on_element_change)
            self.elements_buttons[text].pack(side=Tkinter.LEFT)

        #self.button = Tkinter.Button(self.control_frame, text="QUIT", fg="red", command=self.on_close)
        #self.button.pack(side=Tkinter.LEFT)

        self.progress_bar = Meter(self.control_frame)
        self.progress_bar.pack(side=Tkinter.LEFT)
        self.control_frame.pack()

        frame = Tkinter.Frame(self)
        self.stats_scrollbar = Tkinter.Scrollbar(frame, orient=Tkinter.VERTICAL)
        self.stats = Tkinter.Listbox(frame, height=5, yscrollcommand=self.stats_scrollbar.set, font=('courier',10,'normal'), selectmode=Tkinter.EXTENDED)
        self.stats_scrollbar.pack(side=Tkinter.RIGHT, fill=Tkinter.Y)
        self.stats.pack(fill=Tkinter.BOTH, expand=1)
        frame.pack(fill=Tkinter.X)

    def create_statusbar(self):
        ## create a statusbar
        self.status = StatusBar(self)
        self.status.pack(side=Tkinter.BOTTOM, fill=Tkinter.X)
        self.status.set("hi!")


    def draw_figure_for_first_time(self):
        """ Redraws the figure
        """
        self.axes.grid(True, which="both")
        self.axes.set_title("Spectra", fontsize="10")
        self.axes.set_xlabel("wavelength (nm)", fontsize="10")
        self.axes.set_ylabel("flux", fontsize="10")

        for spec in self.spectra:
            if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
            self.active_spectrum = spec
            self.draw_active_spectrum()
        self.draw_regions("continuum")
        self.draw_regions("lines")
        self.draw_regions("segments")
        self.canvas.draw()

    # Draw all elements in regions array and creates widgets in regions widget array
    # Elements cann be "continuum", "lines" or "segments"
    def draw_regions(self, elements):
        for r in self.region_widgets[elements]:
            r.disconnect_and_remove()
        self.region_widgets[elements] = []

        for r in self.regions[elements]:
            if elements == "lines":
                region = CustomizableRegion(self, "lines", r['wave_base'], r['wave_top'], mark_position=r['wave_peak'], note_text=r['note'])
            else:
                region = CustomizableRegion(self, elements, r['wave_base'], r['wave_top'])
            region.connect()
            self.region_widgets[elements].append(region)


    # Check if exists a spectrum with that name, in that case add a suffix
    def get_name(self, name):
        num_repeated = 0
        max_num = 0
        for spec in self.spectra:
            if spec.name.startswith(name):
                try:
                    # Does it has already a suffix?
                    num = int(spec.name.split(self.dupiclated_name_separator)[-1])
                    if num > max_num:
                        max_num = num
                except ValueError as e:
                    pass
                num_repeated += 1

        if num_repeated > 0: # There are repeated names
            if len(name.split(self.dupiclated_name_separator)) == 0:
                name = name + self.dupiclated_name_separator + str(max_num+1) # Add identificator number + 1
            else:
                name = name.split(self.dupiclated_name_separator)[0] + self.dupiclated_name_separator + str(max_num+1) # Add identificator number + 1
        return name

    def get_color(self):
        # Look for a free color
        free_color_found = None
        if len(self.spectra) < len(self.spectra_colors):
            for color in self.spectra_colors:
                discard_color = False
                for spec in self.spectra:
                    if spec.color == color:
                        discard_color = True
                        break
                if discard_color:
                    continue
                else:
                    free_color_found = color
                    break

        if free_color_found == None:
            good_color = False
            while not good_color:
                # Random number converted to hexadecimal
                random_num = np.random.randint(0, 16777215)
                random_color = "#%x" % random_num
                # Test that the generated color is correct
                try:
                    rgba = matplotlib.colors.colorConverter.to_rgba(random_color)
                    good_color = True
                except ValueError as e:
                    pass

            free_color_found = random_color

        return free_color_found


    def on_license(self):
        license = """Spectra Visual Editor is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

Spectra Visual Editor is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with Spectra Visual Editor.  If not, see:

www.gnu.org/licenses/"""
        self.info("SVE License", license)

    def on_about(self):
        description = """Spectra Visual Editor is a tool for the treatment of spectrum files in order to identify lines, continuum regions and determine radial velocities among other options.
"""
        if "generate_spectrum" in dir(sve):
            description += """
The generation of synthetic spectrum is done thanks to:

SPECTRUM a Stellar Spectral Synthesis Program
(C) Richard O. Gray 1992 - 2010 Version 2.76e
"""
        self.info("About SVE", description)



    def on_action_change(self):
        self.enable_elements()
        self.stats.delete(0, Tkinter.END) # Delete all

        self.action = self.action_entry.get()
        if self.action in ["Stats", "Create"]:
            self.elements_buttons['Line marks'].config(state=Tkinter.DISABLED)
            if self.elements_entry.get() == "Line marks":
                self.elements_entry.set("continuum")
                self.elements = "continuum"
                self.elements_buttons["Continuum"].select()

    def enable_elements(self):
        for text in ["Continuum", "Lines", "Line marks", "Segments"]:
            self.elements_buttons[text].config(state=Tkinter.NORMAL)

    def on_element_change(self):
        if self.elements_entry.get() == "Lines":
            self.elements = "lines"
            self.subelements = None
        elif self.elements_entry.get() == "Line marks":
            self.elements = "lines"
            self.subelements = "marks"
        elif self.elements_entry.get() == "Segments":
            self.elements = "segments"
            self.subelements = None
        else:
            self.elements = "continuum"
            self.subelements = None

    def on_close(self):
        if self.operation_in_progress:
            msg = "There is an operation in progress, are you sure you want to exit anyway?"
            title = "Operation in progress"
            if not self.question(title, msg):
                return

        # Check if there is any spectrum not saved
        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"] or spectra_not_saved:
            msg = "Are you sure you want to exit without saving the regions/spectra?"
            title = "Changes not saved"
            if self.question(title, msg):
                if self.samp_manager != None:
                    self.samp_manager.shutdown()
                self.frame.quit() # stops mainloop
                self.frame.destroy()  # this is necessary on Windows to prevent
                                      # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        else:
            if self.samp_manager != None:
                self.samp_manager.shutdown()
            self.frame.quit() # stops mainloop
            self.frame.destroy()  # this is necessary on Windows to prevent
                                  # Fatal Python Error: PyEval_RestoreThread: NULL tstate




    def on_close_all_spectra(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        some_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                some_not_saved = True
                break

        if some_not_saved:
            msg = "Are you sure you want to close ALL the spectra?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            msg = "Are you sure you want to close ALL the spectra?"
            title = "Close all the spectra"
            if not self.question(title, msg):
                return

        for spec in self.spectra:
            self.axes.lines.remove(spec.plot_id)
            # Remove errors if they exists
            if spec != None and spec.errors_plot_id1 != None:
                self.axes.lines.remove(spec.errors_plot_id1)
                spec.errors_plot_id1 = None
            if spec != None and spec.errors_plot_id2 != None:
                self.axes.lines.remove(spec.errors_plot_id2)
                spec.errors_plot_id2 = None
            # Remove fitted continuum if it exists
            if spec != None and spec.continuum_plot_id != None:
                self.axes.lines.remove(spec.continuum_plot_id)
                spec.continuum_plot_id = None
                spec.continuum_model = None
                spec.continuum_data = None
            # Remove fitted lines if they exist
            for region in self.region_widgets["lines"]:
                if region.line_model.has_key(spec):
                    if region.line_plot_id[spec] != None:
                        self.axes.lines.remove(region.line_plot_id[spec])
                    del region.line_plot_id[spec]
                    del region.line_model[spec]
                    del region.line_extra[spec]
        self.spectra = []
        self.active_spectrum = None
        self.active_spectrum_history = []

        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.update_scale()
        self.flash_status_message("All spectra closed.")

    def on_open_continuum(self):
        if self.check_operation_in_progress():
            return
        self.open_file("continuum")

    def on_open_lines(self):
        if self.check_operation_in_progress():
            return
        self.open_file("lines")

    def on_open_segments(self):
        if self.check_operation_in_progress():
            return
        self.open_file("segments")

    def on_open_spectra(self):
        if self.check_operation_in_progress():
            return
        self.open_file("spectra")

    def open_file(self, elements):

        if elements != "spectra" and self.not_saved[elements]:
            msg = 'Are you sure you want to open a new %s file without saving the current regions?' % elements
            title = 'Changes not saved'
            if not self.question(title, msg):
                self.flash_status_message("Discarded.")
                return


        if elements == "spectra":
            if self.active_spectrum != None and self.active_spectrum.path != None:
                filename = self.active_spectrum.path.split('/')[-1]
                filename_length = len(filename)
                dirname = self.active_spectrum.path[:-filename_length]
            elif self.active_spectrum != None:
                filename = self.active_spectrum.name
                dirname = os.getcwd()
            else:
                filename = ""
                dirname = os.getcwd()
        else:
            if self.filenames[elements] != None:
                filename = self.filenames[elements].split('/')[-1]
                filename_length = len(filename)
                dirname = self.filenames[elements][:-filename_length]
            else:
                filename = ""
                dirname = os.getcwd()
        # create a new thread when a button is pressed

        action_ended = False
        while not action_ended:
            answer_ok = False
            if elements == "spectra":
                ftypes = [('All files', '*'), ('FITS', '*.fit *.fits'), ('Plain text', '*.txt')]
                if sys.platform == "darwin":
                    ftypes = [] # Not working properly in MacOSX
                answer = tkFileDialog.askopenfilenames(title="Open %s..." % elements, initialdir=dirname, filetypes=ftypes, defaultextension=".txt")
                if len(np.unique(answer)) > 0:
                    answer_ok = True
                else:
                    answer_ok = False
            else:
                ftypes = [('All files', '*'), ('Plain text', '*.txt')]
                if sys.platform == "darwin":
                    ftypes = [] # Not working properly in MacOSX
                answer = tkFileDialog.askopenfilename(title="Open %s..." % elements, initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
                if isinstance(answer, basestring) and len(answer) > 0:
                    answer_ok = True
                else:
                    answer_ok = False

            if answer_ok:
                try:
                    if elements == "spectra":
                        paths = np.unique(answer)
                        some_does_not_exists = False
                        for i, path in enumerate(paths):
                            if not os.path.exists(path):
                                msg = 'File %s does not exist.' % os.path.basename(path)
                                title = 'File does not exist'
                                self.error(title, msg)
                                some_does_not_exists = True
                                break
                        if some_does_not_exists:
                            continue # Give the oportunity to select another file name

                        # Remove current continuum from plot if exists
                        self.remove_drawn_continuum_spectrum()

                        # Remove current drawn fitted lines if they exist
                        self.remove_drawn_fitted_lines()

                        for path in paths:
                            # Remove "[A]  " from spectrum name (legend) if it exists
                            if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
                            new_spectrum_data = sve.read_spectrum(path)
                            name = self.get_name(path.split('/')[-1]) # If it already exists, add a suffix
                            color = self.get_color()
                            self.active_spectrum = Spectrum(new_spectrum_data, name, path = path, color=color)
                            self.spectra.append(self.active_spectrum)
                            self.active_spectrum_history.append(self.active_spectrum)
                            self.update_menu_active_spectrum()
                            self.draw_active_spectrum()
                        self.update_scale()

                        if len(paths) == 1:
                            self.flash_status_message("Opened file %s" % paths[0])
                        else:
                            self.flash_status_message("Opened %i spectra files" % len(paths))
                    else:
                        path = answer
                        if not os.path.exists(path):
                            msg = 'File %s does not exist.' % os.path.basename(path)
                            title = 'File does not exist'
                            self.error(title, msg)
                            continue # Give the oportunity to select another file name
                        if elements == "continuum":
                            self.regions[elements] = sve.read_continuum_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        elif elements == "lines":
                            self.regions[elements] = sve.read_line_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        else:
                            # 'segments'
                            self.regions[elements] = sve.read_segment_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        self.flash_status_message("Opened file %s" % path)
                    self.filenames[elements] = path
                    self.canvas.draw()
                    action_ended = True
                except Exception as e:
                    msg = 'A file does not have a compatible format.'
                    title = 'File formatincompatible'
                    self.error(title, msg)
                    continue
            else:
                self.flash_status_message("Discarded.")
                action_ended = True



    def on_save_plot(self):
        if self.check_operation_in_progress():
            return
        file_choices = "PNG (*.png)|*.png"

        if self.active_spectrum != None:
            if self.active_spectrum.path != None:
                filename = self.active_spectrum.path.split('/')[-1] + ".png"
                filename_length = len(filename)
                dirname = self.active_spectrum.path[:-filename_length]
            else:
                filename = self.active_spectrum.name.split(self.dupiclated_name_separator)[0] + ".png"
                dirname = os.getcwd()
        else:
            filename = "SVE_plot_image.png"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            ftypes = [('PNG (*.png)', '*.png')]
            if sys.platform == "darwin":
                ftypes = [] # Not working properly in MacOSX
            answer = tkFileDialog.asksaveasfilename(title="Save plot as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".png")
            if isinstance(answer, basestring) and len(answer) > 0:
                answer_ok = True
            else:
                answer_ok = False

            if answer_ok:
                path = answer
                self.canvas.print_figure(path, dpi=self.dpi)
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                action_ended = True


    def on_save_spectrum(self):
        if self.check_operation_in_progress():
            return

        if not self.check_active_spectrum_exists():
            return

        if self.active_spectrum.path != None:
            filename = self.active_spectrum.path.split('/')[-1]
            filename_length = len(filename)
            dirname = self.active_spectrum.path[:-filename_length]
        else:
            filename = self.active_spectrum.name.split(self.dupiclated_name_separator)[0] + ".txt"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            ftypes = [('All files', '*'), ('FITS', '*.fit *.fits'), ('Plain text', '*.txt'), ('Compressed plain text', '*.gz')]
            if sys.platform == "darwin":
                ftypes = [] # Not working properly in MacOSX
            answer = tkFileDialog.asksaveasfilename(title="Save spectrum as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
            if isinstance(answer, basestring) and len(answer) > 0:
                answer_ok = True
            else:
                answer_ok = False

            if answer_ok:
                path = answer
                self.status_message("Saving %s..." % path)
                # Save, compress if the filename ends with ".gz"
                sve.write_spectrum(self.active_spectrum.data, path)
                self.active_spectrum.not_saved = False
                self.update_title()

                # Change name and path
                self.active_spectrum.path = path
                self.active_spectrum.name = path.split('/')[-1]
                self.active_spectrum.plot_id.set_label("[A] " + self.active_spectrum.name)
                self.update_legend()
                self.canvas.draw()

                # Menu active spectrum
                index = int(self.menu_active_spectrum_num.get())
                self.menu_active_spectrum.entryconfig(index, label=self.active_spectrum.name)

                saved = True
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                saved = False
                action_ended = True

        return saved



    ## Save to a file
    # elements can be "continuum", "lines" or "segments"
    def save_regions(self, elements):
        file_choices = "All|*"
        saved = False
        elements = elements.lower()

        if len(self.region_widgets[elements]) == 0:
            msg = "There is no regions to be saved"
            title = 'Empty regions'
            self.error(title, msg)
            return

        if self.filenames[elements] != None:
            filename = self.filenames[elements].split('/')[-1]
            filename_length = len(filename)
            dirname = self.filenames[elements][:-filename_length]
        else:
            filename = elements + ".txt"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            ftypes = [('All files', '*')]
            if sys.platform == "darwin":
                ftypes = [] # Not working properly in MacOSX
            answer = tkFileDialog.asksaveasfilename(title="Save regions as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
            if isinstance(answer, basestring) and len(answer) > 0:
                answer_ok = True
            else:
                answer_ok = False

            if answer_ok:
                path = answer
                # Already asked by tkFileDialog
                #if os.path.exists(path):
                    #msg = 'Are you sure you want to overwrite the file %s?' % os.path.basename(path)
                    #title = 'File already exists'
                    #if not self.question(title, msg):
                        #continue # Give the oportunity to select a new file name
                self.filenames[elements] = path
                output = open(path, "w")
                # Write header
                if elements == "lines":
                    output.write("wave_peak\twave_base\twave_top\tnote\n")
                else:
                    output.write("wave_base\twave_top\n")

                ## Visible regions: maybe they have been modified, removed or created
                for region in self.region_widgets[elements]:
                    # If the widget is not visible, it has been deleted by the user
                    if region.axvspan.get_visible():
                        if elements == "lines":
                            if region.note != None:
                                note = region.note.get_text()
                            else:
                                note = ""
                            output.write("%.4f" % region.get_wave_peak() + "\t" + "%.4f" % region.get_wave_base() + "\t" + "%.4f" % region.get_wave_top() + "\t" + note + "\n")
                        else:
                            output.write("%.4f" % region.get_wave_base() + "\t" + "%.4f" % region.get_wave_top() + "\n")

                output.close()
                self.regions_saved(elements)
                saved = True
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                saved = False
                action_ended = True

        return saved


    def on_save_continuum_regions(self):
        if self.check_operation_in_progress():
            return
        self.save_regions("continuum")

    def on_save_line_regions(self):
        if self.check_operation_in_progress():
            return
        self.save_regions("lines")

    def on_save_segments(self):
        if self.check_operation_in_progress():
            return
        self.save_regions("segments")

    def on_close_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.not_saved:
            msg = "Are you sure you want to close the spectrum without saving it?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.spectra.remove(self.active_spectrum)
        self.axes.lines.remove(self.active_spectrum.plot_id)


        # Remove errors if they exists
        self.remove_drawn_errors_spectrum()
        # Remove fitted continuum if it exists
        self.remove_drawn_continuum_spectrum()
        # Remove fitted lines if they exist
        self.remove_fitted_lines()
        if len(self.spectra) == 0:
            self.active_spectrum = None
            self.active_spectrum_history = []
        else:
            while True:
                try:
                    self.active_spectrum_history.remove(self.active_spectrum)
                except ValueError:
                    break # All references to current active spectra removed from history
            if len(self.active_spectrum_history) > 0:
                # Activate the previous active spectrum
                self.active_spectrum = self.active_spectrum_history[-1]
            else:
                self.active_spectrum = self.spectra[0]

        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.update_scale()
        self.flash_status_message("Spectrum closed.")

    def remove_drawn_errors_spectrum(self):
        if self.active_spectrum != None and self.active_spectrum.errors_plot_id1 != None:
            self.axes.lines.remove(self.active_spectrum.errors_plot_id1)
            self.active_spectrum.errors_plot_id1 = None
        if self.active_spectrum != None and self.active_spectrum.errors_plot_id2 != None:
            self.axes.lines.remove(self.active_spectrum.errors_plot_id2)
            self.active_spectrum.errors_plot_id2 = None


    def on_motion(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        if self.programmed_flash_status != None: return
        if self.operation_in_progress: return
        self.status_message("Cursor on wavelength %.4f" % event.xdata + " and flux %.4f" % event.ydata)

    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        if self.check_operation_in_progress():
            return

        new_halfwidth = 0.05
        # If the left button is clicked when action "create" is active,
        # create a new region of the selected element type
        # NOTE: line marks only can be created from a line region
        if self.action == "Create" and self.subelements != "marks" and event.button == 1 and event.key == None:
            if self.elements == "continuum":
                region = CustomizableRegion(self, "continuum", event.xdata - new_halfwidth, event.xdata + new_halfwidth)
                region.connect()
                self.region_widgets['continuum'].append(region)
            elif self.elements == "lines":
                note_text = self.ask_value('Note for the new line region:', 'Note', '')
                if note_text == None:
                    # User has clicked "Cancel", stop the line creation
                    return
                region = CustomizableRegion(self, "lines", event.xdata - new_halfwidth, event.xdata + new_halfwidth, mark_position=event.xdata, note_text=note_text)
                region.connect()
                self.region_widgets['lines'].append(region)
                self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - new_halfwidth) + " to " + "%.4f" % (event.xdata + new_halfwidth))
            elif self.elements == "segments":
                factor = 100
                region = CustomizableRegion(self, "segments", event.xdata - factor*new_halfwidth, event.xdata + factor*new_halfwidth)
                region.connect()
                self.region_widgets['segments'].append(region)
                self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - factor*new_halfwidth) + " to " + "%.4f" % (event.xdata + factor*new_halfwidth))
            self.canvas.draw()
            self.regions_changed(self.elements)

    # Change saved status and title for the concrete element
    def regions_saved(self, element_type):
        self.not_saved[element_type] = False
        self.update_title()

    ############################################################################
    def show_marks(self):
        elements = "lines"
        for r in self.region_widgets[elements]:
            r.show_mark()

    def hide_marks(self):
        elements = "lines"
        for r in self.region_widgets[elements]:
            r.hide_mark()


    def update_menu_active_spectrum(self):
        # Remove everything from the list (but keep commands, last one is show errors)
        base = self.menu_active_spectrum.index("Show errors in plot") + 2
        self.menu_active_spectrum.delete(base, Tkinter.END)

        if len(self.spectra) == 0:
            # No spectra loaded
            self.menu_active_spectrum.add_radiobutton(label="None", variable=self.menu_active_spectrum_num, value=str(1), indicatoron=0, state=Tkinter.DISABLED)
        else:
            # Add as many options as spectra
            for i in np.arange(len(self.spectra)):
                if self.spectra[i] == None:
                    continue

                self.menu_active_spectrum.add_radiobutton(label=self.spectra[i].name, variable=self.menu_active_spectrum_num, value=str(i+1), indicatoron=1, command=self.on_change_active_spectrum)

                if self.active_spectrum == self.spectra[i]:
                    self.menu_active_spectrum_num.set(str(i+1))

        if len(self.spectra) > 0:
            for menu, label in self.spectrum_function_items:
                index = menu.index(label)
                menu.entryconfig(index, state=Tkinter.NORMAL)
        else:
            for menu, label in self.spectrum_function_items:
                index = menu.index(label)
                menu.entryconfig(index, state=Tkinter.DISABLED)

    def draw_active_spectrum(self):
        if self.active_spectrum != None:
            # Remove spectrum plot if exists
            if self.active_spectrum.plot_id != None:
                self.axes.lines.remove(self.active_spectrum.plot_id)

            # zorder = 1, always in the background
            self.active_spectrum.plot_id = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'], lw=1, color=self.active_spectrum.color, linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor=self.active_spectrum.color, zorder=1, label="[A] "+self.active_spectrum.name)[0]

            # Draw errors
            if self.show_errors.get():
                if self.active_spectrum.errors_plot_id1 != None:
                    self.axes.lines.remove(self.active_spectrum.errors_plot_id1)
                if self.active_spectrum.errors_plot_id2 != None:
                    self.axes.lines.remove(self.active_spectrum.errors_plot_id2)

                # zorder = 1, always in the background
                self.active_spectrum.errors_plot_id1 = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'] + self.active_spectrum.data['err'], lw=1, color=self.active_spectrum.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]
                self.active_spectrum.errors_plot_id2 = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'] - self.active_spectrum.data['err'], lw=1, color=self.active_spectrum.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]

        # Put a legend to the right of the current axis
        self.update_legend()

    def update_legend(self):
        if len(self.spectra) == 0:
            self.axes.legend_ = None
        else:
            leg = self.axes.legend(loc='upper right', bbox_to_anchor=(1.4, 1.1), ncol=1, shadow=False)
            ltext  = leg.get_texts()
            plt.setp(ltext, fontsize='8')


    def add_stats(self, k, v):
        self.stats.insert(Tkinter.END, "%-40s: %s" % (str(k), str(v)))

    ############################################################################
    def on_change_active_spectrum(self):
        if self.check_operation_in_progress():
            return

        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Remove "[A]  " from spectrum name (legend)
        self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
        # Change active spectrum
        i = int(self.menu_active_spectrum_num.get())-1
        self.active_spectrum = self.spectra[i]
        self.active_spectrum_history.append(self.active_spectrum)
        self.draw_active_spectrum()

        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.canvas.draw()
        self.flash_status_message("Active spectrum: %s" % self.active_spectrum.name)


    def on_flash_status_off(self):
        self.status.set('')
        # Progress bar to zero
        self.update_progress(0)
        self.programmed_flash_status = None


    def on_duplicate_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum != None and self.active_spectrum.plot_id != None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
        new_spectrum_data = self.active_spectrum.data.copy()
        name = self.get_name(self.active_spectrum.name) # If it already exists, add a suffix
        path = self.active_spectrum.path
        color = self.get_color()
        self.active_spectrum = Spectrum(new_spectrum_data, name, path = path, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.active_spectrum_history.append(self.active_spectrum)
        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        #self.draw_continuum_spectrum()
        #self.draw_fitted_lines()
        self.update_scale()
        self.flash_status_message("Spectrum duplicated!")


    def on_receive_spectrum(self, new_spectrum_data, name):
        #if self.check_operation_in_progress():
            #return
        #self.after_idle(self.set_operation_in_progress)
        #self.after_idle(self.flash_status_message, "Receiving spectrum from external application...", progress=False)
        self.after_idle(self.on_receive_spectrum_thread, new_spectrum_data, name)

    def on_receive_spectrum_thread(self, new_spectrum_data, name):
        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum != None and self.active_spectrum.plot_id != None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
        name = self.get_name(name) # If it already exists, add a suffix
        color = self.get_color()
        self.active_spectrum = Spectrum(new_spectrum_data, name, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.active_spectrum_history.append(self.active_spectrum)
        self.update_menu_active_spectrum()
        self.draw_active_spectrum()
        self.update_scale()

        self.flash_status_message("Received file %s" % name)
        self.operation_in_progress = False
        return

    def on_send_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        if self.samp_manager == None:
            return

        if not self.samp_manager.is_connected():
            msg = "No compatible external application can be detected because SVE is not connected to any SAMP hub.\n\n* A SAMP hub can be created by using TOPCAT application."
            title = 'Connection not available'
            self.error(title, msg)
            self.flash_status_message("Not connected to any SAMP hub.")
            return

        ids, names, as_tables = self.samp_manager.get_subscribers()

        if len(names) == 0:
            msg = "No compatible external SAMP application has been detected."
            title = 'Send spectrum error'
            self.error(title, msg)
            self.flash_status_message("No compatible SAMP application detected.")
            return

        key = "SendSpectrumDialog"
        if not self.dialog.has_key(key):
            self.dialog[key] = SendSpectrumDialog(self, "Send spectrum to external application", names)
        self.dialog[key].show()

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        selected_application = self.dialog[key].results["Application"]
        self.dialog[key].destroy()

        selected_application_index = None
        for i in np.arange(len(names)):
            if selected_application == names[i]:
                selected_application_index = i
                break

        if selected_application_index == None:
            raise Exception("This should not happen")

        self.operation_in_progress = True
        self.flash_status_message("Sending spectrum to %s..." % selected_application)
        self.update_progress(10)

        thread = threading.Thread(target=self.on_send_spectrum_thread, args=(selected_application, ids[selected_application_index], as_tables[selected_application_index]))
        thread.setDaemon(True)
        thread.start()

    def on_send_spectrum_thread(self, application, application_id, as_table):
        errors = False
        try:
            self.samp_manager.broadcast_spectrum(self.active_spectrum.data, self.active_spectrum.name, application_id, as_table=as_table)
        except Exception:
            errors = True
        self.after_idle(self.on_send_spectrum_finnish, application, errors)

    def on_send_spectrum_finnish(self, application, errors):
        if errors:
            self.flash_status_message("Communication failed with %s." % (application))
        else:
            self.flash_status_message("Spectrum sent to %s." % (application))
        self.operation_in_progress = False


    def remove_drawn_errors_spectra(self):
        for spec in self.spectra:
            if spec.errors_plot_id1 != None:
                self.axes.lines.remove(spec.errors_plot_id1)
                spec.errors_plot_id1 = None
            if spec.errors_plot_id2 != None:
                self.axes.lines.remove(spec.errors_plot_id2)
                spec.errors_plot_id2 = None
        self.canvas.draw()

    def draw_errors_spectra(self):
        for spec in self.spectra:
            # Remove continuum plot if exists
            if spec.errors_plot_id1 != None:
                self.axes.lines.remove(spec.errors_plot_id1)
            if spec.errors_plot_id2 != None:
                self.axes.lines.remove(spec.errors_plot_id2)

            # zorder = 1, always in the background
            spec.errors_plot_id1 = self.axes.plot(spec.data['waveobs'], spec.data['flux'] + spec.data['err'], lw=1, color=spec.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]
            spec.errors_plot_id2 = self.axes.plot(spec.data['waveobs'], spec.data['flux'] - spec.data['err'], lw=1, color=spec.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]
            # More beautiful but slower:
            #x = self.axes.fill_between(spec.data['waveobs'], spec.data['flux'] - spec.data['err'], spec.data['flux'] + spec.data['err'], color='#CCCCCC',zorder=1)

        self.canvas.draw()

    def on_show_errors(self):
        if self.check_operation_in_progress():
            self.show_errors.set(False)
            return
        if self.show_errors.get():
            self.draw_errors_spectra()
        else:
            self.remove_drawn_errors_spectra()


    def __update_numpy_arrays_from_widgets(self, elements):
        total_regions = len(self.region_widgets[elements])
        if elements == "lines":
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                    self.regions[elements]['wave_peak'][i] = region.get_wave_peak()
                    self.regions[elements]['note'][i] = region.get_note_text()
                i += 1
        else:
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_base', float),('wave_top', float)])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                i += 1

    ####################################################################
    def check_operation_in_progress(self):
        if self.operation_in_progress:
            msg = "This action cannot be done while there is another operation in progress"
            title = "Operation in progress"
            self.error(title, msg)
            return True
        return False

    def check_active_spectrum_exists(self):
        if self.active_spectrum == None or len(self.active_spectrum.data['waveobs']) == 0:
            msg = "There is no active spectrum"
            title = "Spectrum not present"
            self.error(title, msg)
            return False
        return True

    def check_continuum_model_exists(self):
        if self.active_spectrum.continuum_model == None:
            msg = "Please, execute a general continuum fit first"
            title = "Continuum model not fitted"
            self.error(title, msg)
            return False
        return True

    def error(self, title, msg):
        tkMessageBox.showerror(title, msg)

    def info(self, title, msg):
        tkMessageBox.showinfo(title, msg)

    def question(self, title, msg):
        answer_yes = tkMessageBox.askquestion(title, msg) == 'yes'
        if not answer_yes:
            self.flash_status_message("Discarded")
        return answer_yes

    def ask_value(self, text, title, default_value):
        response = tkSimpleDialog.askstring(title, text, initialvalue=default_value)
        if isinstance(response, Tkinter.NoneType):
            response = None
        return response

    def update_progress(self, value):
        self.progress_bar.set((1.*value)/100)

    def update_title(self):
        title = "SVE"

        spectra_not_saved = False
        for spec in self.spectra:
            if spec.not_saved:
                spectra_not_saved = True
                break
        if self.not_saved["continuum"] or self.not_saved["lines"] or self.not_saved["segments"] or spectra_not_saved:
            title += "("
            if self.not_saved["continuum"]:
                title += "*continuum"
            if self.not_saved["lines"]:
                title += "*lines"
            if self.not_saved["segments"]:
                title += "*segments"
            for spec in self.spectra:
                if spec.not_saved:
                    title += "*" + spec.name
            title += "* not saved)"
        self.wm_title(title)

    def update_scale(self):
        # If there are line marks, hide them or they will affect the autoscale
        self.hide_marks()
        # Autoscale
        self.axes.relim()
        #self.axes.autoscale_view(tight=None, scalex=False, scaley=True)
        # Save current view
        ylim = self.axes.get_ylim()
        xlim = self.axes.get_xlim()
        # Change view to fit all the ploted spectra
        self.axes.autoscale(enable=True, axis='both', tight=None)
        self.toolbar.update() # Reset history and consider this view as the default one when 'home' is pressed
        self.toolbar.push_current() # Save the view in the history
        # Recover the line marks in case they existed before
        self.show_marks()
        # Recover initial view
      	self.axes.set_ylim(ylim)
      	self.axes.set_xlim(xlim)
        self.toolbar.push_current() # Save the view in the history
        self.canvas.draw()

    # Change saved status and title for the concrete element
    # - Called from CustomizableRegion
    def regions_changed(self, element_type):
        self.not_saved[element_type] = True
        self.update_title()

    # - Called from CustomizableRegion
    def update_stats(self, region):
        if not self.check_active_spectrum_exists():
            return

        self.stats.delete(0, Tkinter.END) # Delete all

        wave_base = region.get_wave_base()
        wave_top = region.get_wave_top()

        wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
        spectrum_window = self.active_spectrum.data[wave_filter]
        num_points = len(spectrum_window['flux'])

        self.add_stats("Number of measures", num_points)

        self.add_stats("Wavelength min.", "%.4f" % wave_base)
        self.add_stats("Wavelength max.", "%.4f" % wave_top)
        wave_diff = wave_top - wave_base
        self.add_stats("Wavelength difference", "%.4f" % wave_diff)

        if region.element_type == "lines":
            self.add_stats("Wavelength peak", "%.4f" % region.get_wave_peak())
            note = region.get_note_text()
            if note != "":
                self.add_stats("Note", note)

        if num_points > 0:
            self.add_stats("Flux min.", "%.6f" % np.min(spectrum_window['flux']))
            self.add_stats("Flux max.", "%.6f" % np.max(spectrum_window['flux']))
            self.add_stats("Flux mean", "%.6f" % np.mean(spectrum_window['flux']))
            self.add_stats("Flux median", "%.6f" % np.median(spectrum_window['flux']))
            self.add_stats("Flux standard deviation", "%.6f" % np.std(spectrum_window['flux']))

        if region.element_type == "lines" and region.line_plot_id.has_key(self.active_spectrum) and region.line_model[self.active_spectrum] != None:
            self.add_stats("Gaussian mean (mu)", "%.4f" % region.line_model[self.active_spectrum].mu())
            self.add_stats("Gaussian amplitude (A)", "%.4f" % region.line_model[self.active_spectrum].A())
            self.add_stats("Gaussian standard deviation (sigma)", "%.4f" % region.line_model[self.active_spectrum].sig())
            self.add_stats("Gaussian base level (mean continuum)", "%.4f" % region.line_model[self.active_spectrum].baseline())
            rms = region.line_model[self.active_spectrum].rms
            self.add_stats("Gaussian fit root mean squeare (RMS)", "%.4f" % rms)

            #from_x = region.get_wave_peak() - 3*regions['sig'][i]
            #to_x = region.get_wave_peak() + 3*regions['sig'][i]
                #regions['integrated_flux'][i] = -1 * line_model.integrate(from_x, to_x)
            wave_base = region.get_wave_base()
            wave_top = region.get_wave_top()
            integrated_flux = -1 * region.line_model[self.active_spectrum].integrate(wave_base, wave_top)
            ew = integrated_flux / region.line_model[self.active_spectrum].baseline()
            self.add_stats("Gaussian fit Equivalent Width (EW)", "%.4f" % ew)

        if region.element_type == "lines" and region.line_extra.has_key(self.active_spectrum) and region.line_extra[self.active_spectrum] != None:
            # Extras (all in string format separated by ;)
            VALD_wave_peak, species, lower_state, upper_state, loggf, fudge_factor, transition_type, rad, stark, waals, ew, element, telluric_wave_peak, telluric_depth = region.line_extra[self.active_spectrum].split(";")
            self.add_stats("VALD element", element)
            self.add_stats("VALD line wavelength", VALD_wave_peak)
            self.add_stats("VALD lower state (cm^-1)", lower_state)
            self.add_stats("VALD upper state (cm^-1)", upper_state)
            self.add_stats("VALD log(gf)", loggf)
            self.add_stats("VALD radiative damping constant", rad)
            self.add_stats("VALD stark damping constant", stark)
            self.add_stats("VALD van der Waals damping constant", waals)
            if telluric_wave_peak != "" and float(telluric_wave_peak) != 0:
                self.add_stats("Tellurics: possibly affected by line at (nm)", telluric_wave_peak)
                self.add_stats("Tellurics: typical line depth", "%.4f" % float(telluric_depth))
            else:
                self.add_stats("Tellurics:", "Unaffected")

        if self.active_spectrum.continuum_model != None:
            if num_points > 0:
                mean_continuum = np.mean(self.active_spectrum.continuum_model(spectrum_window['waveobs']))
                self.add_stats("Continuum mean for the region", "%.4f" % mean_continuum)
            residuals = np.abs(self.active_spectrum.continuum_model.residuals())
            rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
            self.add_stats("Continuum fit root mean square (RMS)", "%.4f" % rms)



    def flash_status_message(self, msg, flash_len_ms=3000, progress=True):
        if self.programmed_flash_status != None:
            self.after_cancel(self.programmed_flash_status)
            self.programmed_flash_status = None

        self.status.set(msg)
        if progress:
            self.update_progress(100)

        self.programmed_flash_status = self.after(flash_len_ms, self.on_flash_status_off)

    def status_message(self, msg):
        if self.programmed_flash_status != None:
            self.after_cancel(self.programmed_flash_status)
            self.programmed_flash_status = None
            # Progress bar to zero
            self.update_progress(0)

        self.status.set(msg)

    def remove_drawn_continuum_spectrum(self):
        if self.active_spectrum != None and self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
            self.active_spectrum.continuum_plot_id = None
            self.active_spectrum.continuum_model = None
            #self.active_spectrum.continuum_data = None
            self.canvas.draw()

    def remove_drawn_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_plot_id.has_key(self.active_spectrum) and region.line_plot_id[self.active_spectrum] != None:
                self.axes.lines.remove(region.line_plot_id[self.active_spectrum])
                region.line_plot_id[self.active_spectrum] = None
        self.canvas.draw()

    def remove_fitted_lines(self):
        self.remove_drawn_fitted_lines()
        # Remove fitted lines if they exist
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum):
                del region.line_plot_id[self.active_spectrum]
                del region.line_model[self.active_spectrum]
                del region.line_extra[self.active_spectrum]
        if self.active_spectrum != None:
            self.active_spectrum.linemasks = None
            self.active_spectrum.abundances = None

    def draw_continuum_spectrum(self):
        if self.active_spectrum == None:
            return

        # Remove continuum plot if exists
        if self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)

        # zorder = 1, always in the background
        if self.active_spectrum.continuum_data != None:
            self.active_spectrum.continuum_plot_id = self.axes.plot(self.active_spectrum.continuum_data['waveobs'], self.active_spectrum.continuum_data['flux'], lw=1, color='green', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]

    def draw_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum) and region.line_model[self.active_spectrum] != None:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()
                wave_filter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
                spectrum_window = self.active_spectrum.data[wave_filter]

                # Get fluxes from model
                line_fluxes = region.line_model[self.active_spectrum](spectrum_window['waveobs'])

                # zorder = 4, above the line region
                line_plot_id = self.axes.plot(spectrum_window['waveobs'], line_fluxes, lw=1, color='red', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=4)[0]
                region.line_plot_id[self.active_spectrum] = line_plot_id
                self.canvas.draw()


    def on_fit_continuum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return


        key = "FitContinuumDialog"
        # Initial recommendation: 1 knot every 10 nm
        nknots = np.max([1, int((np.max(self.active_spectrum.data['waveobs']) - np.min(self.active_spectrum.data['waveobs'])) / 10)])
        if not self.active_spectrum.dialog.has_key(key):
            median_wave_range=0.1
            max_wave_range=1
            self.active_spectrum.dialog[key] = FitContinuumDialog(self, "Properties for fitting continuum", nknots, median_wave_range, max_wave_range)
        self.active_spectrum.dialog[key].show(suggested_nknots=nknots)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        nknots = self.active_spectrum.dialog[key].results["Number of splines"]
        median_wave_range = self.active_spectrum.dialog[key].results["Wavelength step for median selection"]
        max_wave_range = self.active_spectrum.dialog[key].results["Wavelength step for max selection"]
        in_continuum = self.active_spectrum.dialog[key].results["Fit using"] == "Only continuum regions"
        self.active_spectrum.dialog[key].destroy()

        if nknots == None or median_wave_range < 0 or max_wave_range < 0:
            self.flash_status_message("Bad value.")
            return

        if in_continuum and len(self.region_widgets["continuum"]) == 0:
            self.flash_status_message("No continuum regions found.")
            return

        self.operation_in_progress = True
        self.status_message("Fitting continuum...")
        self.update_progress(10)
        thread = threading.Thread(target=self.on_fit_continuum_thread, args=(nknots,), kwargs={'in_continuum':in_continuum, 'median_wave_range':median_wave_range, 'max_wave_range':max_wave_range})
        thread.setDaemon(True)
        thread.start()

    def __get_spectrum_from_model(self, model, spectrum_wave_grid):
        """
        Calculate flux for a wavelength grid using a given model (i.e. continuum fitted model).
        """
        total_points = len(spectrum_wave_grid)
        spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        spectrum['waveobs'] = spectrum_wave_grid
        spectrum['flux'] = model(spectrum['waveobs'])
        spectrum['err'] = np.zeros(total_points)
        return spectrum

    def on_fit_continuum_thread(self, nknots, in_continuum=False, median_wave_range=0.1, max_wave_range=1):
        try:
            if in_continuum:
                self.__update_numpy_arrays_from_widgets("continuum")
                self.active_spectrum.continuum_model = sve.fit_continuum(self.active_spectrum.data, segments=self.regions["continuum"] , nknots=nknots, median_wave_range=median_wave_range, max_wave_range=max_wave_range)
            else:
                self.active_spectrum.continuum_model = sve.fit_continuum(self.active_spectrum.data, nknots=nknots, median_wave_range=median_wave_range, max_wave_range=max_wave_range)
            self.active_spectrum.continuum_data = self.__get_spectrum_from_model(self.active_spectrum.continuum_model, self.active_spectrum.data['waveobs'])

            self.after_idle(self.on_fit_continuum_finish, nknots)
        except Exception as e:
            self.operation_in_progress = False
            self.after_idle(self.flash_status_message, "Not enough data points to fit, reduce the number of nknots or increase the spectrum regions.")


    def on_fit_continuum_finish(self, nknots):
        self.draw_continuum_spectrum()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Continuum fitted with %s knots uniform Spline model." % str(nknots))

    def on_determine_velocity_atomic(self):
        self.on_determine_velocity(relative_to_atomic_data = True, relative_to_telluric_data = False, relative_to_template=False)

    def on_determine_velocity_telluric(self):
        self.on_determine_velocity(relative_to_atomic_data = False, relative_to_telluric_data = True, relative_to_template=False)

    def on_determine_velocity_template(self):
        self.on_determine_velocity(relative_to_atomic_data = False, relative_to_telluric_data = False, relative_to_template=True)

    def on_determine_velocity(self, relative_to_atomic_data = True, relative_to_telluric_data = False, relative_to_template=False, show_previous_results=True):
        if self.check_operation_in_progress():
            return

        # Default values
        if relative_to_atomic_data:
            title = "Velocity profile relative to atomic lines"
            linelist = self.linelist_atomic
            velocity_lower_limit = self.velocity_atomic_lower_limit
            velocity_upper_limit = self.velocity_atomic_upper_limit
            velocity_step = self.velocity_atomic_step
            templates = []
        elif relative_to_telluric_data:
            title = "Velocity profile relative to telluric lines"
            linelist = self.linelist_telluric
            velocity_lower_limit = self.velocity_telluric_lower_limit
            velocity_upper_limit = self.velocity_telluric_upper_limit
            velocity_step = self.velocity_telluric_step
            templates = []
        elif relative_to_template:
            title = "Velocity profile relative to template"
            linelist = self.linelist_atomic
            velocity_lower_limit = self.velocity_template_lower_limit
            velocity_upper_limit = self.velocity_template_upper_limit
            velocity_step = self.velocity_template_step
            templates = ["[Internal template]"]
            # Add as many options as spectra
            for i in np.arange(len(self.spectra)):
                if self.spectra[i] == None:
                    continue
                templates.append(self.spectra[i].name)
        else:
            raise Exception("Velocity should be determined relative to something!")

        key = "VelocityProfileDialog:"+str(relative_to_atomic_data)+str(relative_to_telluric_data)+str(relative_to_template)
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = VelocityProfileDialog(self, title, rv_lower_limit=velocity_lower_limit, rv_upper_limit=velocity_upper_limit, rv_step=velocity_step, templates=templates)
            self.active_spectrum.dialog[key].show()
        elif show_previous_results:
            if relative_to_template:
                self.active_spectrum.dialog[key].show(updated_templates=templates)
            else:
                self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        rv_lower_limit = self.active_spectrum.dialog[key].results["Velocity lower limit (km/s)"]
        rv_upper_limit = self.active_spectrum.dialog[key].results["Velocity upper limit (km/s)"]
        rv_step = self.active_spectrum.dialog[key].results["Velocity steps (km/s)"]
        if relative_to_template:
            template = self.active_spectrum.dialog[key].results["Cross-correlate with"]
        else:
            template = None
        self.active_spectrum.dialog[key].destroy()

        if rv_lower_limit == None or rv_upper_limit == None or rv_step == None:
            self.flash_status_message("Bad value.")
            return

        if rv_lower_limit >= rv_upper_limit:
            msg = "Upper velocity limit should be greater than lower limit"
            title = "Velocity error"
            self.error(title, msg)
            return

        if (np.abs(rv_lower_limit) + np.abs(rv_upper_limit)) <= 4*rv_step:
            msg = "Velocity step too big for the established limits"
            title = "Velocity error"
            self.error(title, msg)
            return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_determine_velocity_thread, args=(relative_to_atomic_data, relative_to_telluric_data, relative_to_template, rv_lower_limit, rv_upper_limit, rv_step, template))
        thread.setDaemon(True)
        thread.start()

    def __filter_telluric_lines(self, linelist_telluric, spectrum, velocity_lower_limit, velocity_upper_limit):
        """
        Select telluric lines inside a given velocity limits (km/s).
        """
        # Light speed in vacuum
        c = 299792458.0 # m/s

        ## Select telluric lines of interest
        # Limit to region of interest
        wmin = spectrum['waveobs'][0]
        wmax = spectrum['waveobs'][-1]
        delta_wmin = wmin * (velocity_lower_limit / (c/1000.0))
        delta_wmax = wmax * (velocity_upper_limit / (c/1000.0))
        wfilter = (linelist_telluric['wave_peak'] <= wmax + delta_wmax) & (linelist_telluric['wave_peak'] >= wmin + delta_wmin)
        linelist = linelist_telluric[wfilter]
        # Discard not fitted lines
        rfilter = linelist['rms'] == 9999
        linelist = linelist[~rfilter]
        # Discard too deep or too small lines
        rfilter = (linelist['depth'] <= 0.9) & (linelist['depth'] >= 0.01)
        linelist = linelist[rfilter]
        # Discard outliers FWHM in km/s (which is not wavelength dependent)
        telluric_fwhm = (c / (linelist['wave_peak'] / linelist['fwhm'])) / 1000.0 # km/s
        fwhm_selected, fwhm_selected_filter = sve.sigma_clipping(telluric_fwhm, meanfunc=np.median)
        linelist = linelist[fwhm_selected_filter]
        return linelist


    def on_determine_velocity_thread(self, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, velocity_lower_limit, velocity_upper_limit, velocity_step, template):
        self.after_idle(self.status_message, "Determining velocity...")

        if relative_to_atomic_data:
            if self.linelist_atomic == None:
                vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
                self.linelist_atomic = sve.read_VALD_linelist(vald_linelist_file, minimum_depth=0.0)
            linelist = self.linelist_atomic
        elif relative_to_telluric_data:
            if self.linelist_telluric == None:
                telluric_lines_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
                self.linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)
            linelist = self.__filter_telluric_lines(self.linelist_telluric, self.active_spectrum.data, velocity_lower_limit, velocity_upper_limit)
        elif relative_to_template:
            linelist = None
        else:
            raise Exception("Velocity should be determined relative to something!")

        if relative_to_template:
            if template == "[Internal template]":
                # Internal template (solar type)
                if self.active_spectrum.velocity_profile_internal_template == None:
                    self.active_spectrum.velocity_profile_internal_template = sve.read_spectrum(resource_path("input/spectra/synthetic/synth_5777.0_4.44_0.02_2.0.txt.gz"))
                    if self.linelist_telluric == None:
                        telluric_lines_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
                        self.linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)

                    # - Filter regions that may be affected by telluric lines
                    #   (only the 25% of the deepest ones)
                    dfilter = self.linelist_telluric['depth'] > np.percentile(self.linelist_telluric['depth'], 75)
                    tfilter = sve.create_filter_for_regions_affected_by_tellurics(self.active_spectrum.velocity_profile_internal_template['waveobs'], \
                                                self.linelist_telluric[dfilter], min_velocity=-30.0, max_velocity=30.0)
                    self.active_spectrum.velocity_profile_internal_template['flux'][tfilter] = 0.0
                template_spectrum = self.active_spectrum.velocity_profile_internal_template
            else:
                # Search template to be used by its name
                for i in np.arange(len(self.spectra)):
                    if self.spectra[i] == None:
                        continue
                    if self.spectra[i].name == template:
                        template_spectrum = self.spectra[i].data
                        break
            xcoord, fluxes, errors = sve.build_velocity_profile(self.active_spectrum.data, template=template_spectrum, lower_velocity_limit=velocity_lower_limit, upper_velocity_limit=velocity_upper_limit, velocity_step=velocity_step, frame=self)
        else:
            xcoord, fluxes, errors = sve.build_velocity_profile(self.active_spectrum.data, linelist=linelist, lower_velocity_limit=velocity_lower_limit, upper_velocity_limit=velocity_upper_limit, velocity_step=velocity_step, frame=self)

        self.after_idle(self.on_determine_velocity_finish, xcoord, fluxes, errors, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, linelist)

    def on_determine_velocity_finish(self, xcoord, fluxes, errors, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, linelist):
        # Modelize
        if relative_to_atomic_data or relative_to_template:
            models = sve.modelize_velocity_profile(xcoord, fluxes, errors)
            if len(models) > 1:
                accept = sve.select_good_velocity_profile_models(models, xcoord, fluxes)
                if len(models[accept]) == 0:
                    models = models[:1]
                else:
                    models = models[accept]
        else:
            models = sve.modelize_velocity_profile(xcoord, fluxes, errors, only_one_peak=True)

        if len(models) == 0:
            fwhm = 0.0
            telluric_fwhm = 0.0
            R = 0.0
            velocity = 0.0
            self.flash_status_message("Velocity could not be determined!")
        else:
            # Resolving power
            c = 299792458.0 # m/s
            fwhm = models[0].fwhm()[0] # km/s (because xcoord is already velocity)
            if relative_to_atomic_data or relative_to_template:
                telluric_fwhm = 0.0
                R = np.int(c/(1000.0*fwhm))
            else:
                # If telluric lines have been used, we can substract its natural FWHM
                # so that we get the real resolution of the instrument (based on the difference in FWHM)
                c = 299792458.0 # m/s
                telluric_fwhm = np.mean((c / (linelist['wave_peak'] / linelist['fwhm'])) / 1000.0) # km/s
                R = np.int(c/(1000.0*np.round(fwhm - telluric_fwhm, 2)))
            # Velocity
            velocity = np.round(models[0].mu(), 2) # km/s
            if relative_to_atomic_data:
                self.active_spectrum.velocity_atomic = velocity
                self.active_spectrum.resolution_atomic = R
            elif relative_to_telluric_data:
                self.active_spectrum.velocity_telluric = velocity
                self.active_spectrum.resolution_telluric = R
            elif relative_to_template:
                self.active_spectrum.velocity_template = velocity
            # A positive velocity indicates the distance between the objects is or was increasing;
            # A negative velocity indicates the distance between the source and observer is or was decreasing.
            self.flash_status_message("Velocity determined: " + str(velocity) + " km/s")
        self.operation_in_progress = False

        key = "VelocityProfileDialog:"+str(relative_to_atomic_data)+str(relative_to_telluric_data)+str(relative_to_template)
        self.active_spectrum.dialog[key].register(xcoord, fluxes, errors, models, telluric_fwhm=telluric_fwhm)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return
        else:
            # Recalculate
            self.active_spectrum.dialog[key].destroy()
            self.on_determine_velocity(relative_to_atomic_data=relative_to_atomic_data, relative_to_telluric_data = relative_to_telluric_data, relative_to_template = relative_to_template, show_previous_results=False)



    def on_fit_lines(self):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return

        total_regions = len(self.region_widgets["lines"])
        if total_regions == 0:
            msg = 'There are no line regions to fit.'
            title = "No line regions"
            self.error(title, msg)
            return

        key = "FitLinesDialog"
        vel_telluric = self.active_spectrum.velocity_telluric
        vel_atomic = self.active_spectrum.velocity_atomic
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = FitLinesDialog(self, "Fit lines and cross-match with VALD linelist", vel_telluric, vel_atomic)
        self.active_spectrum.dialog[key].show(updated_vel_telluric=vel_telluric, updated_vel_atomic=vel_atomic)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        vel_atomic = self.active_spectrum.dialog[key].results["Velocity respect to atomic lines (km/s)"]
        vel_telluric = self.active_spectrum.dialog[key].results["Velocity respect to telluric lines (km/s)"]
        self.active_spectrum.dialog[key].destroy()

        if vel_atomic == None or vel_telluric == None:
            self.flash_status_message("Bad value.")
            return

        self.active_spectrum.velocity_atomic = vel_atomic
        self.active_spectrum.velocity_telluric = vel_telluric

        # Remove drawd lines if they exist
        self.remove_drawn_fitted_lines()


        self.operation_in_progress = True
        self.status_message("Fitting lines...")
        thread = threading.Thread(target=self.on_fit_lines_thread, args=(vel_atomic, vel_telluric))
        thread.setDaemon(True)
        thread.start()

    def on_fit_lines_thread(self, vel_atomic, vel_telluric):
        self.__update_numpy_arrays_from_widgets("lines")
        vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")
        telluric_linelist_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
        linemasks = sve.fit_lines(self.regions["lines"], self.active_spectrum.data, self.active_spectrum.continuum_model, vel_atomic, vel_telluric, vald_linelist_file, chemical_elements_file, molecules_file, telluric_linelist_file, discard_gaussian=False, discard_voigt=True, smoothed_spectrum=self.active_spectrum.data, frame=self)
        # Exclude lines that have not been successfully cross matched with the atomic data
        # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
        rejected_by_atomic_line_not_found = (linemasks['VALD_wave_peak'] == 0)
        linemasks = linemasks[~rejected_by_atomic_line_not_found]

        self.after_idle(self.on_fit_lines_finnish, linemasks)

    def on_fit_lines_finnish(self, linemasks, conserve_previous_regions=True):
        elements = "lines"

        diff_num_regions = len(self.region_widgets["lines"]) - len(linemasks)
        if conserve_previous_regions and diff_num_regions > 0:
            # Find regions that has been discarded due to a bad fit or other reason
            # in order to recover them
            recovered_regions = np.recarray((diff_num_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            i = 0
            for region in self.region_widgets["lines"]:
                lost_region = True
                for line in linemasks:
                    if region.get_wave_base() == line['wave_base'] and region.get_wave_top() == line['wave_top']:
                        lost_region = False
                        continue
                if lost_region:
                    recovered_regions[i]['wave_base'] = region.get_wave_base()
                    recovered_regions[i]['wave_top'] = region.get_wave_top()
                    recovered_regions[i]['wave_peak'] = region.get_wave_peak()
                    recovered_regions[i]['note'] = region.get_note_text()
                    i += 1

        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements, check_not_saved=False)

        self.active_spectrum.linemasks = linemasks

        if linemasks != None and len(linemasks) > 0:
            total_regions = len(linemasks)
            line_regions = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            line_regions['wave_peak'] = linemasks['mu']
            # If no edge limit improvement has been applied (i.e. fit_lines does not do it)
            line_regions['wave_base'] = linemasks['wave_base']
            line_regions['wave_top'] = linemasks['wave_top']
            line_regions['note'] = linemasks['element']

            # Save the data in the note, separated by commas (only the first element will be shown in the GUI)
            line_models = []
            line_extras = []
            i = 0
            for line in linemasks:
                line_extra = str(line['VALD_wave_peak']) + ";" + str(line['species']) + ";"
                line_extra = line_extra + str(line['lower state (cm^-1)']) + ";" + str(line['upper state (cm^-1)']) + ";"
                line_extra = line_extra + str(line['log(gf)']) + ";" + str(line['fudge factor']) + ";"
                line_extra = line_extra + str(line['transition type']) + ";" + str(line['rad']) + ";"
                line_extra = line_extra + str(line['stark']) + ";" + str(line['waals']) + ";"
                line_extra = line_extra + str(line['ew']) + ";" + line['element'] + ";"
                line_extra = line_extra + str(line['telluric_wave_peak']) + ";" + str(line['telluric_depth'])
                line_extras.append(line_extra)
                line_model = sve.GaussianModel(baseline=line['baseline'], A=line['A'], sig=line['sig'], mu=line['mu'])
                line_model.rms = line['rms']
                line_models.append(line_model)
                if line["telluric_wave_peak"] != 0:
                    line_regions['note'][i] += "*"
                i += 1

            if conserve_previous_regions and diff_num_regions > 0:
                line_regions = np.hstack((line_regions, recovered_regions))

            self.regions[elements] = line_regions
            self.draw_regions(elements)

            # Fitted Gaussian lines
            # * The not fitted but recovered lines will be avoided by the "diff_num_regions + i"
            for i in np.arange(len(line_models)):
                self.region_widgets["lines"][i].line_model[self.active_spectrum] = line_models[i]
                self.region_widgets["lines"][i].line_extra[self.active_spectrum] = line_extras[i]
            self.draw_fitted_lines()

            self.not_saved[elements] = True
            self.update_title()
            self.flash_status_message("%i line masks found/fitted!" % (len(line_regions)))
        else:
            self.flash_status_message("No line masks found/fitted!")
        self.canvas.draw()
        self.operation_in_progress = False


    def remove_regions(self, elements, check_not_saved=True):
        if check_not_saved and self.not_saved[elements]:
            msg = "Are you sure you want to remove this regions without saving them?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        if elements == "continuum":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        elif elements == "lines":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
        else:
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        self.draw_regions(elements)
        self.not_saved[elements] = False
        self.update_title()
        self.filenames[elements] = None
        self.canvas.draw()


    def remove_continuum_spectrum(self):
        self.remove_drawn_continuum_spectrum()
        self.active_spectrum.continuum_model = None
        self.active_spectrum.continuum_data = None

    def on_remove_fitted_continuum(self):
        if self.check_operation_in_progress():
            return
        self.remove_continuum_spectrum()

    def on_remove_fitted_lines(self):
        if self.check_operation_in_progress():
            return
        self.remove_fitted_lines()

    def remove_regions(self, elements, check_not_saved=True):
        if check_not_saved and self.not_saved[elements]:
            msg = "Are you sure you want to remove this regions without saving them?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        if elements == "continuum":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        elif elements == "lines":
            self.regions[elements] = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
        else:
            self.regions[elements] = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
        self.draw_regions(elements)
        self.not_saved[elements] = False
        self.update_title()
        self.filenames[elements] = None
        self.canvas.draw()

    def on_remove_continuum_regions(self):
        if self.check_operation_in_progress():
            return
        elements = "continuum"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned continuum regions")


    def on_remove_line_masks(self):
        if self.check_operation_in_progress():
            return
        elements = "lines"
        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements)
        self.flash_status_message("Cleaned line masks")

    def on_remove_segments(self):
        if self.check_operation_in_progress():
            return
        elements = "segments"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned segments")



    def on_find_continuum(self):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return

        if self.not_saved["continuum"]:
            msg = "Are you sure you want to find new continuum regions without saving the current ones?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        key = "FindContinuumDialog"
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = FindContinuumDialog(self, "Properties for finding continuum regions", self.find_continuum_regions_wave_step, self.find_continuum_regions_sigma, self.find_continuum_regions_max_continuum_diff)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        resolution = None
        fixed_wave_step = self.active_spectrum.dialog[key].results["Check for regions of minimum size"]
        sigma = self.active_spectrum.dialog[key].results["Maximum standard deviation"]
        max_continuum_diff = self.active_spectrum.dialog[key].results["Maximum fitted continuum difference (%)"]
        in_segments = self.active_spectrum.dialog[key].results["Look for continuum regions in"] == "Only inside segments"
        self.active_spectrum.dialog[key].destroy()

        if fixed_wave_step == None or sigma == None or max_continuum_diff == None:
            self.after_idel(self.flash_status_message, "Bad value.")
            return
        # Save values
        self.find_continuum_regions_wave_step = fixed_wave_step
        self.find_continuum_regions_sigma = sigma
        self.find_continuum_regions_max_continuum_diff = max_continuum_diff
        # Convert from % to over 1
        max_continuum_diff = max_continuum_diff / 100

        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            self.after_idle(self.flash_status_message, "No segments found.")
            return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_find_continuum_thread, args=(resolution, sigma, max_continuum_diff), kwargs={'fixed_wave_step':fixed_wave_step, 'in_segments':in_segments})
        thread.setDaemon(True)
        thread.start()

    def update_numpy_arrays_from_widgets(self, elements):
        total_regions = len(self.region_widgets[elements])
        if elements == "lines":
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                    self.regions[elements]['wave_peak'][i] = region.get_wave_peak()
                    self.regions[elements]['note'][i] = region.get_note_text()
                i += 1
        else:
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_base', float),('wave_top', float)])
            i = 0
            for region in self.region_widgets[elements]:
                if region.axvspan.get_visible():
                    self.regions[elements]['wave_base'][i] = region.get_wave_base()
                    self.regions[elements]['wave_top'][i] = region.get_wave_top()
                i += 1

    def on_find_continuum_thread(self, resolution, sigma, max_continuum_diff, fixed_wave_step=None, in_segments=False):
        self.after_idle(self.status_message, "Finding continuum regions...")
        if in_segments:
            self.update_numpy_arrays_from_widgets("segments")
            continuum_regions = sve.find_continuum(self.active_spectrum.data, resolution, segments=self.regions["segments"], max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)
        else:
            continuum_regions = sve.find_continuum(self.active_spectrum.data, resolution, max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)

        self.after_idle(self.on_find_continuum_finish, continuum_regions)

    def on_find_continuum_finish(self, continuum_regions):
        elements = "continuum"
        self.regions[elements] = continuum_regions
        self.draw_regions(elements)
        self.not_saved[elements] = True
        self.update_title()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Automatic finding of continuum regions ended.")

    def on_find_lines(self):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return

        if self.not_saved["lines"]:
            msg = "Are you sure you want to find new line masks without saving the current ones?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        if self.active_spectrum.resolution_telluric == 0.0:
            R = self.active_spectrum.resolution_atomic
        else:
            R = self.active_spectrum.resolution_telluric

        vel_atomic = self.active_spectrum.velocity_atomic
        vel_telluric = self.active_spectrum.velocity_telluric


        key = "FindLinesDialog"
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = FindLinesDialog(self, "Properties for finding line masks", self.find_lines_min_depth, self.find_lines_max_depth, vel_atomic=vel_atomic, vel_telluric=vel_telluric, resolution=R, elements="Fe 1, Fe 2")
        self.active_spectrum.dialog[key].show(updated_vel_atomic=vel_atomic, updated_vel_telluric=vel_telluric)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        min_depth = self.active_spectrum.dialog[key].results["Minimum depth (% of the continuum)"]
        max_depth = self.active_spectrum.dialog[key].results["Maximum depth (% of the continuum)"]
        elements = self.active_spectrum.dialog[key].results["Select elements (comma separated)"]
        resolution = self.active_spectrum.dialog[key].results["Resolution"]
        vel_atomic = self.active_spectrum.dialog[key].results["Velocity respect to atomic lines (km/s)"]
        vel_telluric = self.active_spectrum.dialog[key].results["Velocity respect to telluric lines (km/s)"]
        discard_tellurics = self.active_spectrum.dialog[key].results["Discard affected by tellurics"] == 1
        in_segments = self.active_spectrum.dialog[key].results["Look for line masks in"] == "Only inside segments"
        self.active_spectrum.dialog[key].destroy()

        if max_depth == None or min_depth == None or resolution == None or vel_atomic == None or vel_telluric == None or max_depth <= min_depth or max_depth <= 0 or min_depth < 0 or resolution <= 0:
            self.after_idle(self.flash_status_message, "Bad value.")
            return
        ## Save values
        self.find_lines_max_depth = max_depth
        self.find_lines_min_depth = min_depth
        self.active_spectrum.velocity_atomic = vel_atomic
        self.active_spectrum.velocity_telluric = vel_telluric
        self.active_spectrum.resolution_telluric = resolution
        self.active_spectrum.resolution_atomic = resolution

        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            self.after_idle(self.flash_status_message, "No segments found.")
            return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_find_lines_thread, args=(max_depth, min_depth, elements, resolution, vel_atomic, vel_telluric, discard_tellurics), kwargs={'in_segments':in_segments})
        thread.setDaemon(True)
        thread.start()

    def on_find_lines_thread(self, max_depth, min_depth, elements, resolution, vel_atomic, vel_telluric, discard_tellurics, in_segments=False):
        if in_segments:
            # Select spectrum from regions
            self.update_numpy_arrays_from_widgets("segments")
            spectrum = None
            for region in self.region_widgets["segments"]:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()

                wfilter = np.logical_and(self.active_spectrum.data['waveobs'] >= wave_base, self.active_spectrum.data['waveobs'] <= wave_top)
                if spectrum == None:
                    spectrum = self.active_spectrum.data[wfilter]
                else:
                    spectrum = np.hstack((spectrum, self.active_spectrum.data[wfilter]))
        else:
            spectrum = self.active_spectrum.data

        logging.info("Smoothing spectrum...")
        self.after_idle(self.status_message, "Smoothing spectrum...")
        smoothed_spectrum = sve.convolve_spectrum(spectrum, 2*resolution, frame=self)


        self.after_idle(self.status_message, "Generating line masks, fitting gaussians and matching VALD lines...")
        logging.info("Generating line masks, fitting gaussians and matching VALD lines...")
        vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")
        telluric_linelist_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
        linemasks = sve.find_linemasks(spectrum, self.active_spectrum.continuum_model, vald_linelist_file, chemical_elements_file, molecules_file, telluric_linelist_file, minimum_depth=min_depth, maximum_depth=max_depth, smoothed_spectrum=smoothed_spectrum, discard_gaussian = False, discard_voigt = True, vel_atomic=vel_atomic, vel_telluric=vel_telluric, frame=self)

        # If no peaks found, just finnish
        if linemasks == None or len(linemasks) == 0:
            self.after_idle(self.on_find_lines_finish, None)
            return

        logging.info("Applying filters to discard bad line masks...")
        self.after_idle(self.status_message, "Applying filters to discard bad line masks...")

        rejected_by_atomic_line_not_found = (linemasks['VALD_wave_peak'] == 0)
        rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)

        discarded = linemasks['wave_peak'] <= 0 # All to false

        # In case it is specified, select only given elements
        if elements != "":
            elements = elements.split(",")
            select = linemasks['element'] == "NONE" # All to false
            for element in elements:
                select = np.logical_or(select, linemasks['element'] == element.strip().capitalize())
            discarded = np.logical_or(discarded, np.logical_not(select))

        if discard_tellurics:
            discarded = np.logical_or(discarded, rejected_by_telluric_line)

        # Exclude lines that have not been successfully cross matched with the atomic data
        # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
        discarded = np.logical_or(discarded, rejected_by_atomic_line_not_found)

        if in_segments:
            # Identify linemasks with too big wavelength range (probably star in between segments)
            wave_diff = (linemasks['wave_top'] - linemasks['wave_base']) / (linemasks['top'] - linemasks['base'])
            #wave_diff_selected, wave_diff_selected_filter = sve.sigma_clipping(wave_diff[~discarded], sig=3, meanfunc=np.median) # Discard outliers
            wave_diff_selected, wave_diff_selected_filter = sve.interquartile_range_filtering(wave_diff[~discarded], k=1.5) # Discard outliers
            accepted_index = np.arange(len(wave_diff))[~discarded][wave_diff_selected_filter]
            rejected_by_wave_gaps = wave_diff < 0 # Create an array of booleans
            rejected_by_wave_gaps[:] = True # Initialize
            rejected_by_wave_gaps[accepted_index] = False
            discarded = np.logical_or(discarded, rejected_by_wave_gaps)

        linemasks = linemasks[~discarded]
        total_regions = len(linemasks)

        # If no regions found, just finnish
        if total_regions == 0:
            self.after_idle(self.on_find_lines_finish, None, None, None)
            return


        self.after_idle(self.on_find_lines_finish, linemasks)

    def on_find_lines_finish(self, linemasks):
        self.on_fit_lines_finnish(linemasks, conserve_previous_regions=False)



    def on_correct_velocity_atomic(self):
        vel_type = "atomic lines"
        default_vel = self.active_spectrum.velocity_atomic
        self.on_correct_vel(vel_type, default_vel)

    def on_correct_velocity_telluric(self):
        vel_type = "telluric lines"
        default_vel = self.active_spectrum.velocity_telluric
        self.on_correct_vel(vel_type, default_vel)

    def on_correct_velocity_template(self):
        vel_type = "template"
        default_vel = self.active_spectrum.velocity_template
        self.on_correct_vel(vel_type, default_vel)

    def on_correct_vel(self, vel_type, default_vel):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        key = "CorrectVelocityDialog:"+vel_type
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = CorrectVelocityDialog(self, "Velocity correction", vel_type, default_vel)
        self.active_spectrum.dialog[key].show(updated_vel=default_vel)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return


        velocity = self.active_spectrum.dialog[key].results["Velocity relative to %s (km/s)" % vel_type]
        in_regions = self.active_spectrum.dialog[key].results["Apply correction on"] == "Regions"
        self.active_spectrum.dialog[key].destroy()

        if velocity == None:
            self.flash_status_message("Bad value.")
            return


        self.status_message("Correcting " + vel_type + " velocity...")

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        if in_regions:
            for elements in ["lines", "continuum", "segments"]:
                self.update_numpy_arrays_from_widgets(elements)
                if len(self.regions[elements]) > 0:
                    if elements == "lines":
                        self.regions[elements] = sve.correct_velocity_regions(self.regions[elements], velocity, with_peak=True)
                    else:
                        self.regions[elements] = sve.correct_velocity_regions(self.regions[elements], velocity)
                    self.draw_regions(elements)
                    self.not_saved[elements] = False
        else:
            # Modify velocities according to the correction
            self.active_spectrum.velocity_atomic -= velocity
            self.active_spectrum.velocity_telluric -= velocity
            self.active_spectrum.velocity_template -= velocity
            # Correct
            self.active_spectrum.data = sve.correct_velocity(self.active_spectrum.data, velocity)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a " + vel_type + " velocity correction of %s." % velocity)


    def on_determine_barycentric_vel(self):
        key = "DetermineBarycentricCorrectionDialog"
        if not self.dialog.has_key(key):
            self.dialog[key] = DetermineBarycentricCorrectionDialog(self,  "Barycentric velocity determination", "15/02/2012", "00:00:00", "19:50:46.99", "08:52:5.96")
        self.dialog[key].show()

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        date_string = self.dialog[key].results["Date (DD/MM/YYY)"]
        time_string = self.dialog[key].results["Time (HH:MM:SS)"]
        ra = self.dialog[key].results["Right ascension (HH:MM:SS)"]
        dec = self.dialog[key].results["Declination (DD:MM:SS)"]
        self.dialog[key].destroy()

        try:
            day, month, year = map(float, date_string.split("/"))
            hours, minutes, seconds = map(float, time_string.split(":"))
            ra_hours, ra_minutes, ra_seconds = map(float, ra.split(":"))
            dec_degrees, dec_minutes, dec_seconds = map(float, dec.split(":"))
        except:
            msg = 'Some input values are not in the expected data format.'
            title = "Bad values"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        if None in [day, month, year, hours, minutes, seconds, ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds]:
            self.flash_status_message("Bad value.")
            return

        # Project velocity toward star
        self.barycentric_vel = sve.calculate_barycentric_velocity_correction((year, month, day, hours, minutes, seconds), (ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds))

        msg = "Barycentric velocity determined: " + str(self.barycentric_vel) + " km/s"
        title = "Barycentric velocity"
        self.info(title, msg)
        self.flash_status_message(msg)

    def on_estimate_snr(self):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.snr != None:
            msg = "Previous estimated SNR: %.2f. Re-estimate again? " % self.active_spectrum.snr
            title = "Signal-to-Noise Ratio"
            if not self.question(title, msg):
                return

        key = "EstimateSNRDialog"
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = EstimateSNRDialog(self, "Properties for estimating SNR", num_points=10)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        num_points = self.active_spectrum.dialog[key].results["Number of points"]
        wave_step = self.active_spectrum.dialog[key].results["Wavelength step (resampling)"]
        estimate_from_flux = self.active_spectrum.dialog[key].results["Estimate SNR"] == "Directly from reported errors"
        self.active_spectrum.dialog[key].destroy()

        if num_points == None or wave_step == None:
            self.flash_status_message("Bad value.")
            return

        if estimate_from_flux:
            self.operation_in_progress = True
            thread = threading.Thread(target=self.on_estimate_snr_thread, args=(num_points, wave_step))
            thread.setDaemon(True)
            thread.start()
        else:
            efilter = self.active_spectrum.data['err'] > 0
            spec = self.active_spectrum.data[efilter]
            if len(spec) > 1:
                estimated_snr = np.median(spec['flux'] / spec['err'])
                self.on_estimate_snr_finnish(estimated_snr)
            else:
                msg = 'All value errors are set to zero or negative numbers'
                title = "SNR estimation error"
                self.error(title, msg)


    def on_estimate_snr_thread(self, num_points, wave_step):
        self.after_idle(self.status_message, "Resampling spectrum...")
        self.after_idle(self.update_progress, 10)
        xaxis = np.arange(np.min(self.active_spectrum.data["waveobs"]), np.max(self.active_spectrum.data["waveobs"]), wave_step)
        resampled_spectrum_data = sve.resample_spectrum(self.active_spectrum.data, xaxis, frame=self)

        self.after_idle(self.status_message, "Estimating SNR for the whole spectrum...")
        estimated_snr = sve.estimate_snr(self.active_spectrum.data['flux'], num_points=num_points, frame=self)
        self.after_idle(self.on_estimate_snr_finnish, estimated_snr)

    def on_estimate_snr_finnish(self, estimated_snr):
        self.active_spectrum.snr = estimated_snr
        msg = "Estimated SNR: %.2f" % self.active_spectrum.snr
        title = "Signal-to-Noise Ratio"
        self.info(title, msg)
        self.flash_status_message(msg)
        self.operation_in_progress = False

    def on_estimate_errors(self):
        if self.check_operation_in_progress():
            return

        if np.all(self.active_spectrum.data['err'] > 0.0) or np.all(self.active_spectrum.data['err'] < np.max(self.active_spectrum.data['flux'])):
            mean_err = np.mean(self.active_spectrum.data['err'])
            max_err = np.max(self.active_spectrum.data['err'])
            min_err = np.min(self.active_spectrum.data['err'])
            msg = "Spectrum seems to already have errors (mean: %.4f, max: %.4f, min: %.4f). Are you sure to overwrite them with new estimates? " % (mean_err, max_err, min_err)
            title = "Spectrum's errors"
            if not self.question(title, msg):
                return

        key = "EstimateErrorsDialog"
        if not self.active_spectrum.dialog.has_key(key):
            snr = 10.0
            self.active_spectrum.dialog[key] = EstimateErrorsDialog(self, "Estimate spectrum errors", snr)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        snr = self.active_spectrum.dialog[key].results["SNR (Signal-to-Noise Ratio)"]
        self.active_spectrum.dialog[key].destroy()

        if snr == None or snr <= 0:
            msg = "SNR should be greater than zero"
            title = "SNR error"
            self.error(title, msg)
            return

        self.active_spectrum.data['err'] = self.active_spectrum.data['flux'] / snr
        self.active_spectrum.not_saved = True

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        msg = "Spectrum's errors estimated!"
        self.flash_status_message(msg)


    def on_degrade_resolution(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        key = "DegradeResolutionDialog"
        # Give priority to the resolution determined by the velocity profile relative to telluric lines
        if self.active_spectrum.resolution_telluric == 0.0:
            R = self.active_spectrum.resolution_atomic
        else:
            R = self.active_spectrum.resolution_telluric
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = DegradeResolutionDialog(self, "Degrade spectrum resolution", R, R/2)
        self.active_spectrum.dialog[key].show(updated_from_resolution=R, updated_to_resolution=R/2)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        from_resolution = self.active_spectrum.dialog[key].results["Initial resolution"]
        to_resolution = self.active_spectrum.dialog[key].results["Final resolution"]
        self.active_spectrum.dialog[key].destroy()

        if from_resolution == None or to_resolution == None or (from_resolution != 0 and from_resolution <= to_resolution) or from_resolution < 0 or to_resolution <= 0:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to degrade its resolution now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_degrade_resolution_thread, args=(from_resolution, to_resolution,))
        thread.setDaemon(True)
        thread.start()

    def on_degrade_resolution_thread(self, from_resolution, to_resolution):
        self.after_idle(self.status_message, "Degrading spectrum resolution...")
        if from_resolution == 0:
            # Smooth
            convolved_spectrum = sve.convolve_spectrum(self.active_spectrum.data, to_resolution, frame=self)
        else:
            convolved_spectrum = sve.convolve_spectrum(self.active_spectrum.data, from_resolution, to_resolution=to_resolution, frame=self)
        self.after_idle(self.on_degrade_resolution_finnish, convolved_spectrum, from_resolution, to_resolution)

    def on_degrade_resolution_finnish(self, convolved_spectrum, from_resolution, to_resolution):
        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        self.active_spectrum.data = convolved_spectrum
        self.active_spectrum.not_saved = True
        self.active_spectrum.resolution_telluric = to_resolution
        self.active_spectrum.resolution_atomic = to_resolution

        self.draw_active_spectrum()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Spectrum degreaded from resolution %.2f to %.2f" % (from_resolution, to_resolution))
        self.operation_in_progress = False

    def on_clean_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        key = "CleanSpectrumDialog"
        if not self.dialog.has_key(key):
            self.dialog[key] = CleanSpectrumDialog(self, "Clean fluxes and errors", 0.0, 1.2, 0.0, 1.2)
        self.dialog[key].show()

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        filter_by_flux = self.dialog[key].results["Filter by flux"] == 1
        flux_base = self.dialog[key].results["Base flux"]
        flux_top = self.dialog[key].results["Top flux"]
        filter_by_error = self.dialog[key].results["Filter by error"] == 1
        err_base = self.dialog[key].results["Base error"]
        err_top = self.dialog[key].results["Top error"]
        self.dialog[key].destroy()

        if flux_base == None or flux_top == None or flux_top <= flux_base or err_base == None or err_top == None or err_top <= err_base:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to clean it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Cleaning spectrum...")

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        if filter_by_flux and filter_by_error:
            ffilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
            efilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
            wfilter = np.logical_and(ffilter, efilter)
        elif filter_by_flux:
            wfilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
        elif filter_by_error:
            wfilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
        else:
            wfilter = np.logical_not(np.isnan(self.active_spectrum.data['err']))

        if len(self.active_spectrum.data[wfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong flux/error ranges"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        self.active_spectrum.data = self.active_spectrum.data[wfilter]
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum cleaned.")

    def on_clean_tellurics(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        key = "CleanTelluricsDialog"
        rv = self.active_spectrum.velocity_atomic
        if not self.active_spectrum.dialog.has_key(key):
            self.active_spectrum.dialog[key] = CleanTelluricsDialog(self, "Clean telluric regions", rv, -30.0, 30.0, 0.02)
        self.active_spectrum.dialog[key].show(updated_vel=rv)

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return

        rv = self.active_spectrum.dialog[key].results["Radial velocity"]
        min_vel = self.active_spectrum.dialog[key].results["Minimum velocity"]
        max_vel = self.active_spectrum.dialog[key].results["Maximum velocity"]
        min_depth = self.active_spectrum.dialog[key].results["Minimum tellurics depth"]
        self.active_spectrum.dialog[key].destroy()

        if rv == None or min_depth == None or min_vel == None or max_vel == None or min_vel >= max_vel:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to clean it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Cleaning spectrum...")

        self.active_spectrum.velocity_atomic = rv

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        if self.linelist_telluric == None:
            telluric_lines_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
            self.linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)

        # - Filter regions that may be affected by telluric lines
        #dfilter = self.linelist_telluric['depth'] > np.percentile(self.linelist_telluric['depth'], 75) # (only the 25% of the deepest ones)
        dfilter = self.linelist_telluric['depth'] > min_depth
        tfilter = sve.create_filter_for_regions_affected_by_tellurics(self.active_spectrum.data['waveobs'], \
                                    self.linelist_telluric[dfilter], min_velocity=-rv+min_vel, max_velocity=-rv+max_vel)

        if len(self.active_spectrum.data[tfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong ranges"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return
        #self.active_spectrum.data['flux'][tfilter] = 0.0
        #self.active_spectrum.data['err'][tfilter] = 0.0
        self.active_spectrum.data = self.active_spectrum.data[~tfilter]

        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum cleaned.")

    def on_cut_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        key = "CutSpectrumDialog"
        if not self.dialog.has_key(key):
            self.dialog[key] = CutSpectrumDialog(self, "Wavelength range reduction", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2))
        self.dialog[key].show()

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        wave_base = self.dialog[key].results["Base wavelength"]
        wave_top = self.dialog[key].results["Top wavelength"]
        in_segments = self.dialog[key].results["Consider"] == "Segments"
        self.dialog[key].destroy()

        if not in_segments and (wave_base == None or wave_top == None or wave_top <= wave_base):
            self.flash_status_message("Bad value.")
            return

        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            self.after_idle(self.flash_status_message, "No segments found.")
            return

        if not in_segments:
            wfilter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
        else:
            # Build wavelength points from regions
            wfilter = None
            for region in self.region_widgets["segments"]:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()

                if wfilter == None:
                    wfilter = np.logical_and(self.active_spectrum.data['waveobs'] >= wave_base, self.active_spectrum.data['waveobs'] <= wave_top)
                else:
                    wfilter = np.logical_or(wfilter, np.logical_and(self.active_spectrum.data['waveobs'] >= wave_base, self.active_spectrum.data['waveobs'] <= wave_top))

        if len(self.active_spectrum.data[wfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong wavelength range"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to cut it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Cutting spectrum...")

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        self.active_spectrum.data = self.active_spectrum.data[wfilter]
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum cutted.")

    def on_resample_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        steps = self.active_spectrum.data['waveobs'][1:] - self.active_spectrum.data['waveobs'][:-1]
        # Round to the first decimal found
        median_step = np.median(steps)
        median_step = np.round(median_step, decimals=np.max([0, int(-1*np.floor(np.log10(median_step)))]))
        mean_step = np.mean(steps)
        mean_step = np.round(mean_step, decimals=np.max([0, int(-1*np.floor(np.log10(mean_step)))]))
        min_step = np.min(steps)
        min_step = np.round(min_step, decimals=np.max([0, int(-1*np.floor(np.log10(min_step)))]))
        max_step = np.max(steps)
        max_step = np.round(max_step, decimals=np.max([0, int(-1*np.floor(np.log10(max_step)))]))

        key = "ResampleSpectrumDialog"
        if not self.dialog.has_key(key):
            self.dialog[key] = ResampleSpectrumDialog(self, "Resample spectrum", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2), 0.001, median_step, mean_step, min_step, max_step)
        self.dialog[key].show(updated_median_step=median_step, updated_mean_step=mean_step, updated_min_step=min_step, updated_max_step=max_step)

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        wave_base = self.dialog[key].results["Base wavelength"]
        wave_top = self.dialog[key].results["Top wavelength"]
        wave_step = self.dialog[key].results["Wavelength step"]
        method = self.dialog[key].results["Method"]
        self.dialog[key].destroy()

        if wave_base == None or wave_top == None or wave_top <= wave_base or wave_step <= 0:
            self.flash_status_message("Bad value.")
            return
        # It is not necessary to check if base and top are out of the current spectrum,
        # the resample_spectrum function can deal it

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to proceed with the resampling anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.operation_in_progress = True

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        thread = threading.Thread(target=self.on_resample_spectrum_thread, args=(wave_base, wave_top, wave_step, method,))
        thread.setDaemon(True)
        thread.start()

    def on_resample_spectrum_thread(self, wave_base, wave_top, wave_step, method):
        # Homogenize
        self.after_idle(self.status_message, "Resampling spectrum...")
        self.after_idle(self.update_progress, 10)
        xaxis = np.arange(wave_base, wave_top, wave_step)
        resampled_spectrum_data = sve.resample_spectrum(self.active_spectrum.data, xaxis, method=method, frame=self)
        self.active_spectrum.data = resampled_spectrum_data
        self.after_idle(self.on_resample_spectrum_finnish)

    def on_resample_spectrum_finnish(self):
        self.active_spectrum.not_saved = True

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum resampled.")
        self.operation_in_progress = False

    def on_combine_spectra(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        wave_base = None
        wave_top = None
        for spec in self.spectra:
            new_top = np.max(spec.data['waveobs'])
            new_base = np.min(spec.data['waveobs'])
            if wave_base == None:
                wave_base = new_base
                wave_top = new_top
                continue
            if wave_base > new_base:
                wave_base = new_base
            if wave_top < new_top:
                wave_top = new_top

        # This dialog is not needed to be saved
        dialog = CombineSpectraDialog(self, "Resample & combine spectrum", np.round(wave_base, 2), np.round(wave_top, 2), 0.001)
        dialog.show()

        if dialog.results == None:
            # Cancel
            dialog.destroy()
            return

        wave_base = dialog.results["Base wavelength"]
        wave_top = dialog.results["Top wavelength"]
        wave_step = dialog.results["Wavelength step"]
        operation = dialog.results["Operation"]
        dialog.destroy()

        operation_median = False
        operation_mean = False
        operation_subtract = False
        operation_add = False
        operation_divide = False
        if operation == "Median":
            operation_median = True
        elif operation == "Mean":
            operation_mean = True
        elif operation == "Subtract":
            operation_subtract = True
        elif operation == "Add":
            operation_add = True
        elif operation == "Divide":
            operation_divide = True
        else:
            raise Exception("Unknown combine operation")


        if wave_base == None or wave_top == None or wave_top <= wave_base or wave_step <= 0:
            self.flash_status_message("Bad value.")
            return
        # It is not necessary to check if base and top are out of the current spectrum,
        # the resample_spectrum function can deal it

        # Check if there is any spectrum not saved
        #spectra_not_saved = False
        #for spec in self.spectra:
            #if spec.not_saved:
                #spectra_not_saved = True
                #break
        #if spectra_not_saved:
            #msg = "Some of the spectra have not been saved, are you sure you want to combine them all without saving previously?"
            #title = "Changes not saved"
            #if not self.question(title, msg):
                #return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_combine_spectra_thread, args=(wave_base, wave_top, wave_step, operation_median, operation_mean, operation_subtract, operation_add, operation_divide))
        thread.setDaemon(True)
        thread.start()

    def on_combine_spectra_thread(self, wave_base, wave_top, wave_step, operation_median, operation_mean, operation_subtract, operation_add, operation_divide):
        # Homogenize
        i = 0
        total = len(self.spectra)
        xaxis = np.arange(wave_base, wave_top, wave_step)
        resampled_spectra = []
        active = 0
        for spec in self.spectra:
            self.after_idle(self.status_message, "Resampling spectrum %i of %i (%s)..." % (i+1, total, spec.name))
            if spec == self.active_spectrum:
                active = i
            resampled_spectrum_data = sve.resample_spectrum(spec.data, xaxis, frame=self)
            resampled_spectra.append(resampled_spectrum_data)
            i += 1

        # Combine
        self.after_idle(self.status_message, "Combining spectra...")
        self.after_idle(self.update_progress, 10)

        total_spectra = len(resampled_spectra)
        total_wavelengths = len(resampled_spectra[0]['waveobs'])
        if operation_median or operation_mean:
            matrix = np.empty((total_spectra, total_wavelengths))

            # Create a matrix for an easy combination
            i = 0
            for spec in resampled_spectra:
                matrix[i] = spec['flux']
                i += 1

            # Standard deviation
            std = np.std(matrix, axis=0)

            if operation_mean:
                # Mean fluxes
                combined_spectrum = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
                combined_spectrum['waveobs'] = xaxis
                combined_spectrum['flux'] = np.mean(matrix, axis=0)
                combined_spectrum['err'] = std
                combined_spectrum_name = "Cumulative_mean_spectrum"
            else:
                # Median fluxes
                combined_spectrum = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
                combined_spectrum['waveobs'] = xaxis
                combined_spectrum['flux'] = np.median(matrix, axis=0)
                combined_spectrum['err'] = std
                combined_spectrum_name = "Cumulative_median_spectrum"
        elif operation_subtract:
            flux = np.zeros(total_wavelengths)
            err = np.zeros(total_wavelengths)
            i = 0
            for spec in resampled_spectra:
                if i == active:
                    flux = flux + spec['flux']
                else:
                    flux = flux - spec['flux']
                # Error propagation assuming that they are independent
                err = np.sqrt(np.power(err,2) + np.power(spec['err'],2))
                i += 1
            combined_spectrum = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
            combined_spectrum['waveobs'] = xaxis
            combined_spectrum['flux'] = flux
            combined_spectrum['err'] = err
            combined_spectrum_name = "Subtracted_spectrum"
        elif operation_add:
            flux = np.zeros(total_wavelengths)
            err = np.zeros(total_wavelengths)
            for spec in resampled_spectra:
                flux = flux + spec['flux']
                # Error propagation assuming that they are independent
                err = np.sqrt(np.power(err,2) + np.power(spec['err'],2))
            combined_spectrum = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
            combined_spectrum['waveobs'] = xaxis
            combined_spectrum['flux'] = flux
            combined_spectrum['err'] = err
            combined_spectrum_name = "Added_spectrum"
        elif operation_divide:
            flux = np.zeros(total_wavelengths)
            err = np.zeros(total_wavelengths)
            i = 0
            for spec in resampled_spectra:
                logging.warn("Division by zero may occur")
                if i == 0:
                    err = spec['err'].copy()
                    flux = spec['flux'].copy()
                else:
                    # Error propagation assuming that they are independent
                    err = np.sqrt(np.power(flux / err, 2) + np.power(spec['flux'] / spec['err'], 2)) * 0.7
                    if i == active:
                        flux = flux * spec['flux']
                    else:
                        flux = flux * (1. / spec['flux'])
                i += 1
            combined_spectrum = np.recarray((total_wavelengths, ), dtype=[('waveobs', float),('flux', float),('err', float)])
            combined_spectrum['waveobs'] = xaxis
            combined_spectrum['flux'] = flux
            combined_spectrum['err'] = err
            combined_spectrum_name = "Divided_spectrum"

        # Free memory
        for spec in resampled_spectra:
            del spec
        del resampled_spectra

        self.after_idle(self.on_combine_spectra_finnish, combined_spectrum, combined_spectrum_name)

    def on_combine_spectra_finnish(self, combined_spectrum, combined_spectrum_name):
        # Remove "[A]  " from spectrum name (legend)
        self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        # Add plot
        name = self.get_name(combined_spectrum_name) # If it already exists, add a suffix
        color = self.get_color()
        self.active_spectrum = Spectrum(combined_spectrum, name, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.active_spectrum_history.append(self.active_spectrum)
        self.draw_active_spectrum()

        self.update_menu_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectra combined.")
        self.operation_in_progress = False


    def on_continuum_normalization(self):
        if not self.check_active_spectrum_exists():
            return
        if not self.check_continuum_model_exists():
            return
        if self.check_operation_in_progress():
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to divide all the fluxes by the fitted continuum anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            msg = "This operation is going to divide all the fluxes of the active spectrum by the fitted continuum. Are you sure?"
            title = "Confirmation"
            if not self.question(title, msg):
                return

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        self.status_message("Continuum normalization...")
        self.active_spectrum.data['flux'] /= self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'])
        self.active_spectrum.data['err'] /= self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'])
        self.active_spectrum.not_saved = True

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        pass


    def on_operate_spectrum(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        flux = self.active_spectrum.data['flux']
        waveobs = self.active_spectrum.data['waveobs']
        err = self.active_spectrum.data['err']
        self.safe_operations['flux'] = flux
        self.safe_operations['waveobs'] = waveobs
        self.safe_operations['err'] = err

        key = "OperateSpectrumDialog"
        if not self.dialog.has_key(key):
            operation_waveobs = self.operation_waveobs
            operation_flux = self.operation_flux
            operation_err = self.operation_err
            self.dialog[key] = OperateSpectrumDialog(self, "Apply mathematical expression", self.safe_operations_description, operation_waveobs, operation_flux, operation_err)
        self.dialog[key].show()

        if self.dialog[key].results == None:
            # Cancel
            self.dialog[key].destroy()
            return

        operation_waveobs = self.dialog[key].results["Wavelengths ="]
        operation_flux = self.dialog[key].results["Fluxes ="]
        operation_err = self.dialog[key].results["Errors ="]
        self.dialog[key].destroy()

        operate_waveobs = True
        operate_flux = True
        operate_err = True
        if operation_waveobs == "waveobs":
            operate_waveobs = False
        if operation_flux == "flux":
            operate_flux = False
        if operation_err == "err":
            operate_err = False

        if not operate_waveobs and not operate_flux and not operate_err:
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to modify it anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Operating...")

        try:
            if operate_waveobs:
                new_waveobs = eval(operation_waveobs,{"__builtins__":None},self.safe_operations)
                if len(waveobs) != len(new_waveobs):
                    raise Exception("Invalid operation!")
            if operate_flux:
                new_flux = eval(operation_flux,{"__builtins__":None},self.safe_operations)
                if len(waveobs) != len(new_flux):
                    raise Exception("Invalid operation!")
            if operate_err:
                new_err = eval(operation_err,{"__builtins__":None},self.safe_operations)
                if len(waveobs) != len(new_err):
                    raise Exception("Invalid operation!")
        except Exception:
            self.flash_status_message("Invalid operation.")
            return


        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        if operate_waveobs:
            self.active_spectrum.data['waveobs'] = new_waveobs
            self.operation_waveobs = operation_waveobs
        if operate_flux:
            self.active_spectrum.data['flux'] = new_flux
            self.operation_flux = operation_flux
        if operate_err:
            self.active_spectrum.data['err'] = new_err
            self.operation_err = operation_err
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Operation succesfully executed.")

    def on_synthesize(self):
        if self.check_operation_in_progress():
            return
        if "generate_spectrum" in dir(sve):
            if self.active_spectrum != None:
                wave_base = np.round(np.min(self.active_spectrum.data['waveobs']), 2)
                wave_top = np.round(np.max(self.active_spectrum.data['waveobs']), 2)
            else:
                wave_base = 515.0 # Magnesium triplet region
                wave_top = 525.0
                wave_top = 517.0
            teff = 5777.0
            logg = 4.44
            MH = 0.02
            macroturbulence = 0.0
            vsini = 2.0
            limb_darkening_coeff = 0.0
            microturbulence_vel = 2.0
            #resolution = 47000
            resolution = 0
            wave_step = 0.001

            key = "SendSpectrumDialog"
            if not self.dialog.has_key(key):
                self.dialog[key] = SyntheticSpectrumDialog(self, "Synthetic spectrum generator", wave_base, wave_top, wave_step, resolution, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff)
            self.dialog[key].show()

            if self.dialog[key].results == None:
                # Cancel
                self.dialog[key].destroy()
                return

            teff = self.dialog[key].results["Effective temperature (K)"]
            logg = self.dialog[key].results["Surface gravity (log g)"]
            MH = self.dialog[key].results["Metallicity [Fe/H]"]
            microturbulence_vel = self.dialog[key].results["Microturbulence velocity (km/s)"]
            macroturbulence = self.dialog[key].results["Macroturbulence velocity (km/s)"]
            vsini = self.dialog[key].results["Rotation (v sin(i)) (km/s)"]
            limb_darkening_coeff = self.dialog[key].results["Limb darkening coefficient"]
            resolution = self.dialog[key].results["Resolution"]
            wave_base = self.dialog[key].results["Wavelength min (nm)"]
            wave_top = self.dialog[key].results["Wavelength max (nm)"]
            wave_step = self.dialog[key].results["Wavelength step (nm)"]
            selected_atmosphere_models = self.dialog[key].results["Model atmosphere"]
            selected_linelist = self.dialog[key].results["Line list"]
            in_segments = self.dialog[key].results["Generate spectrum for"] == "Segments"
            in_lines = self.dialog[key].results["Generate spectrum for"] == "Line masks"

            self.dialog[key].destroy()

            linelist_file = resource_path("input/linelists/SPECTRUM/" + selected_linelist + "/300_1100nm.lst")
            abundances_file = resource_path("input/abundances/" + selected_atmosphere_models + "/stdatom.dat")

            if in_segments:
                elements_type = "segments"
            if in_lines:
                elements_type = "lines"

            if teff == None or logg == None or MH == None or microturbulence_vel == None or resolution == None or wave_base == None or wave_top == None:
                self.flash_status_message("Bad value.")
                return

            if not self.modeled_layers_pack.has_key(selected_atmosphere_models):
                logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.modeled_layers_pack[selected_atmosphere_models] = sve.load_modeled_layers_pack(resource_path('input/atmospheres/' + selected_atmosphere_models + '/modeled_layers_pack.dump'))

            if not sve.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of theatmospheric models."
                title = 'Out of the atmospheric models'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            # Load SPECTRUM linelist
            if not linelist_file in self.linelist_SPECTRUM.keys():
                self.linelist_SPECTRUM[linelist_file] = sve.read_SPECTRUM_linelist(linelist_file)
            linelist = self.linelist_SPECTRUM[linelist_file]

            # Load SPECTRUM abundances
            if not abundances_file in self.abundances_SPECTRUM.keys():
                self.abundances_SPECTRUM[abundances_file] = sve.read_SPECTRUM_abundances(abundances_file)
            abundances = self.abundances_SPECTRUM[abundances_file]

            # Prepare atmosphere model
            self.status_message("Interpolating atmosphere model...")
            atmosphere_layers = sve.interpolate_atmosphere_layers(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH)

            if wave_base >= wave_top:
                msg = "Bad wavelength range definition, maximum value cannot be lower than minimum value."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return
            waveobs = np.arange(wave_base, wave_top, wave_step)

            if not in_segments and not in_lines:
                lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
                waveobs_mask = np.ones(len(waveobs)) # Compute fluxes for all the wavelengths
            else:
                #in_segments or in_lines
                if len(self.region_widgets[elements_type]) == 0:
                    self.flash_status_message("No segments/line masks present for synthetic spectrum generation.")
                    return

                # Build wavelength points from regions
                wfilter = None
                for region in self.region_widgets[elements_type]:
                    wave_base = region.get_wave_base()
                    wave_top = region.get_wave_top()

                    if wfilter == None:
                        lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
                        wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
                    else:
                        lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave (A)'] >= wave_base*11., linelist['wave (A)'] <= wave_top*10.))
                        wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
                waveobs_mask = np.zeros(len(waveobs))
                waveobs_mask[wfilter] = 1.0 # Compute fluxes only for selected segments
            linelist = linelist[lfilter]


            # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectrum
            if np.min(waveobs) < 300.0 or np.max(waveobs) > 1100.0:
                # luke.300_1000nm.lst
                msg = "Wavelength range is outside line list for spectrum generation."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            total_points = len(waveobs)

            if total_points < 2:
                msg = "Wavelength range too narrow."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            ## TODO: Control wavelength margins, they should not be bigger/lower thant the linelist

            self.operation_in_progress = True
            self.status_message("Synthesizing spectrum...")
            self.update_progress(10)
            thread = threading.Thread(target=self.on_synthesize_thread, args=(waveobs, waveobs_mask, linelist, abundances, atmosphere_layers, teff, logg, MH, microturbulence_vel,  macroturbulence, vsini, limb_darkening_coeff, resolution, ))
            thread.setDaemon(True)
            thread.start()

    def on_synthesize_thread(self, waveobs, waveobs_mask, linelist, abundances, atmosphere_layers, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, resolution):
        total_points = len(waveobs)

        synth_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        synth_spectrum['waveobs'] = waveobs
        synth_spectrum['err'] = 0.0

        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])

        # waveobs is multiplied by 10.0 in order to be converted from nm to armstrongs
        synth_spectrum['flux'] = sve.generate_spectrum(synth_spectrum['waveobs'], waveobs_mask, atmosphere_layers, teff, logg, MH, linelist=linelist, abundances=abundances, fixed_abundances=fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, R=resolution, verbose=1, update_progress_func=self.update_progress)


        synth_spectrum.sort(order='waveobs') # Make sure it is ordered by wavelength

        # Remove atmosphere model temporary file
        self.after_idle(self.on_synthesize_finnish, synth_spectrum, teff, logg, MH, microturbulence_vel)

    def on_synthesize_finnish(self, synth_spectrum, teff, logg, MH, microturbulence_vel):
        self.operation_in_progress = False

        # Check if synthetic generation has failed
        if np.all(synth_spectrum['flux'] == 0):
            self.flash_status_message("The synthetic spectrum generation has failed for those astrophysical parameters!")
            return

        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum != None and self.active_spectrum.plot_id != None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        # Name: If it already exists, add a suffix
        name = self.get_name("Synth_" + str(teff) + "_" + str(logg) + "_"  + str(MH) + "_" + str(microturbulence_vel))
        color = self.get_color()
        self.active_spectrum = Spectrum(synth_spectrum, name, color=color)

        self.spectra.append(self.active_spectrum)
        self.active_spectrum_history.append(self.active_spectrum)
        self.active_spectrum.not_saved = True
        self.update_title()
        self.update_menu_active_spectrum()
        self.draw_active_spectrum()
        self.update_scale()

        self.canvas.draw()

        self.operation_in_progress = False
        self.flash_status_message("Synthetic spectrum generated!")

    def on_determine_abundances(self, show_previous_results=True):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.linemasks == None:
            msg = "Lines should be fitted first."
            title = 'Lines not fitted'
            self.error(title, msg)
            self.flash_status_message("Not previous fitted lines available.")
            return

        key = "AbundancesDialog"
        if not self.active_spectrum.dialog.has_key(key):
            teff = 5777.0
            logg = 4.44
            MH = 0.0
            microturbulence_vel = 2.0
            self.active_spectrum.dialog[key] = AbundancesDialog(self, "Abundances determination", teff, logg, MH, microturbulence_vel)
            self.active_spectrum.dialog[key].show()
        elif show_previous_results:
            self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
        logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
        MH = self.active_spectrum.dialog[key].results["Metallicity [Fe/H]"]
        microturbulence_vel = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
        selected_atmosphere_models = self.active_spectrum.dialog[key].results["Model atmosphere"]
        abundances_file = resource_path("input/abundances/" + selected_atmosphere_models + "/stdatom.dat")
        self.active_spectrum.dialog[key].destroy()

        if teff == None or logg == None or MH == None or microturbulence_vel == None:
            self.flash_status_message("Bad value.")
            return


        if not self.modeled_layers_pack.has_key(selected_atmosphere_models):
            logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.modeled_layers_pack[selected_atmosphere_models] = sve.load_modeled_layers_pack(resource_path('input/atmospheres/' + selected_atmosphere_models + '/modeled_layers_pack.dump'))

        if not sve.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of theatmospheric models."
            title = 'Out of the atmospheric models'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return

        # Load SPECTRUM abundances
        if not abundances_file in self.abundances_SPECTRUM.keys():
            self.abundances_SPECTRUM[abundances_file] = sve.read_SPECTRUM_abundances(abundances_file)
        abundances = self.abundances_SPECTRUM[abundances_file]

        # Prepare atmosphere model
        self.status_message("Interpolating atmosphere model...")
        atmosphere_layers = sve.interpolate_atmosphere_layers(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH)

        self.operation_in_progress = True
        self.status_message("Determining abundances...")
        self.update_progress(10)

        thread = threading.Thread(target=self.on_determine_abundances_thread, args=(atmosphere_layers, teff, logg, MH, abundances, microturbulence_vel,))
        thread.setDaemon(True)
        thread.start()


    def on_determine_abundances_thread(self, atmosphere_layers, teff, logg, MH, abundances, microturbulence_vel):
        linemasks = self.active_spectrum.linemasks.copy()
        linemasks['ew'] = 1000. * 10. * linemasks['ew'] # From nm to mA
        linemasks['VALD_wave_peak'] = 10 * linemasks['VALD_wave_peak'] # From nm to Angstrom

        abundances, normal_abundances, relative_abundances = sve.determine_abundances(atmosphere_layers, teff, logg, MH, linemasks, abundances, microturbulence_vel = 2.0, verbose=1, update_progress_func=self.update_progress)

        self.after_idle(self.on_determine_abundances_finnish, abundances, normal_abundances, relative_abundances)

    def on_determine_abundances_finnish(self, abundances, normal_abundances, relative_abundances):
        self.flash_status_message("Abundances determined!")
        self.operation_in_progress = False

        self.active_spectrum.abundances = normal_abundances

        key = "AbundancesDialog"
        self.active_spectrum.dialog[key].register(self.active_spectrum.linemasks, abundances)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results == None:
            self.active_spectrum.dialog[key].destroy()
            return
        else:
            # Recalculate
            self.active_spectrum.dialog[key].destroy()
            self.on_determine_abundances(show_previous_results=False)


    def on_determine_parameters(self):
        if self.check_operation_in_progress():
            return
        if "generate_spectrum" in dir(sve):
            teff = 5777.0
            logg = 4.44
            MH = 0.02
            macroturbulence = 0.0
            vsini = 2.0
            limb_darkening_coeff = 0.0
            microturbulence_vel = 2.0
            #resolution = 47000
            resolution = 0
            wave_step = 0.001

            key = "SolverDialog"
            if not self.active_spectrum.dialog.has_key(key):
                self.active_spectrum.dialog[key] = SolverDialog(self, "Determine parameters", resolution, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff)
            self.active_spectrum.dialog[key].show()

            if self.active_spectrum.dialog[key].results == None:
                # Cancel
                self.active_spectrum.dialog[key].destroy()
                return

            teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
            logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
            MH = self.active_spectrum.dialog[key].results["Metallicity [Fe/H]"]
            microturbulence_vel = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
            macroturbulence = self.active_spectrum.dialog[key].results["Macroturbulence velocity (km/s)"]
            vsini = self.active_spectrum.dialog[key].results["Rotation (v sin(i)) (km/s)"]
            limb_darkening_coeff = self.active_spectrum.dialog[key].results["Limb darkening coefficient"]
            resolution = self.active_spectrum.dialog[key].results["Resolution"]
            selected_atmosphere_models = self.active_spectrum.dialog[key].results["Model atmosphere"]
            selected_linelist = self.active_spectrum.dialog[key].results["Line list"]

            free_teff = self.active_spectrum.dialog[key].results["Free Teff"] == 1
            free_logg = self.active_spectrum.dialog[key].results["Free Log(g)"] == 1
            free_MH = self.active_spectrum.dialog[key].results["Free [Fe/H]"] == 1
            free_microturbulence_vel = self.active_spectrum.dialog[key].results["Free Vmic"] == 1
            free_macroturbulence = self.active_spectrum.dialog[key].results["Free Vmac"] == 1
            free_vsini = self.active_spectrum.dialog[key].results["Free vsin(i)"] == 1
            free_limb_darkening_coeff = self.active_spectrum.dialog[key].results["Free limb dark. coeff."] == 1
            free_resolution = self.active_spectrum.dialog[key].results["Free resolution"] == 1

            self.active_spectrum.dialog[key].destroy()

            linelist_file = resource_path("input/linelists/SPECTRUM/" + selected_linelist + "/300_1100nm.lst")
            abundances_file = resource_path("input/abundances/" + selected_atmosphere_models + "/stdatom.dat")

            if teff == None or logg == None or MH == None or microturbulence_vel == None or resolution == None:
                self.flash_status_message("Bad value.")
                return

            if not self.modeled_layers_pack.has_key(selected_atmosphere_models):
                logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.modeled_layers_pack[selected_atmosphere_models] = sve.load_modeled_layers_pack(resource_path('input/atmospheres/' + selected_atmosphere_models + '/modeled_layers_pack.dump'))

            if not sve.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of theatmospheric models."
                title = 'Out of the atmospheric models'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            # Load SPECTRUM linelist
            if not linelist_file in self.linelist_SPECTRUM.keys():
                self.linelist_SPECTRUM[linelist_file] = sve.read_SPECTRUM_linelist(linelist_file)
            linelist = self.linelist_SPECTRUM[linelist_file]

            # Load SPECTRUM abundances
            if not abundances_file in self.abundances_SPECTRUM.keys():
                self.abundances_SPECTRUM[abundances_file] = sve.read_SPECTRUM_abundances(abundances_file)
            abundances = self.abundances_SPECTRUM[abundances_file]

            # Prepare atmosphere model
            self.status_message("Interpolating atmosphere model...")
            atmosphere_layers = sve.interpolate_atmosphere_layers(self.modeled_layers_pack[selected_atmosphere_models], teff, logg, MH)


            # Consider only segments
            elements_type = "segments"
            if len(self.region_widgets[elements_type]) == 0:
                self.flash_status_message("No segments/line masks present for synthetic spectrum generation.")
                return

            # Build wavelength points from regions
            wfilter = None
            for region in self.region_widgets[elements_type]:
                wave_base = region.get_wave_base()
                wave_top = region.get_wave_top()

                if wfilter == None:
                    lfilter = np.logical_and(linelist['wave (A)'] >= wave_base*10., linelist['wave (A)'] <= wave_top*10.)
                    wfilter = np.logical_and(waveobs >= wave_base, waveobs <= wave_top)
                else:
                    lfilter = np.logical_or(lfilter, np.logical_and(linelist['wave (A)'] >= wave_base*11., linelist['wave (A)'] <= wave_top*10.))
                    wfilter = np.logical_or(wfilter, np.logical_and(waveobs >= wave_base, waveobs <= wave_top))
            waveobs_mask = np.zeros(len(waveobs))
            waveobs_mask[wfilter] = 1.0 # Compute fluxes only for selected segments
            linelist = linelist[lfilter]


            # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectrum
            if np.min(waveobs) < 300.0 or np.max(waveobs) > 1100.0:
                # luke.300_1000nm.lst
                msg = "Wavelength range is outside line list for spectrum generation."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            total_points = len(waveobs)

            if total_points < 2:
                msg = "Wavelength range too narrow."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            self.operation_in_progress = True
            self.status_message("Determining parameters...")
            self.update_progress(10)
            thread = threading.Thread(target=self.on_determine_parameters_thread, args=(waveobs, waveobs_mask, linelist, abundances, atmosphere_layers, teff, logg, MH, microturbulence_vel,  macroturbulence, vsini, limb_darkening_coeff, resolution, ))
            thread.setDaemon(True)
            thread.start()

    def on_determine_parameters_thread(self, waveobs, waveobs_mask, linelist, abundances, atmosphere_layers, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, resolution):
        # TODO
        self.after_idle(self.on_determine_parameters_finnish, synth_spectrum, teff, logg, MH, microturbulence_vel)

    def on_determine_parameters_finnish(self, synth_spectrum, teff, logg, MH, microturbulence_vel):
        self.operation_in_progress = False
        # TODO
        self.error("Operation not implemented!")
        self.flash_status_message("Parameters determined!")




