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
import tkinter
import tkinter.messagebox
import tkinter.filedialog
import tkinter.simpledialog
import matplotlib
#matplotlib.use('TkAgg')

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg as FigCanvas, NavigationToolbar2Tk as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt

#import Queue
#from multiprocessing import Queue
from queue import Queue
from queue import Empty
from astropy.io import ascii
import os
import sys

import numpy as np
import threading
import logging


import ispec
from .dialogs import AbundancesDialog
from .dialogs import AdjustLinesDialog
from .dialogs import AddNoiseDialog
from .dialogs import CleanSpectrumDialog
from .dialogs import CleanTelluricsDialog
from .dialogs import CombineSpectraDialog
from .dialogs import CorrectVelocityDialog
from .dialogs import CutSpectrumDialog
from .dialogs import DegradeResolutionDialog
from .dialogs import DetermineBarycentricCorrectionDialog
from .dialogs import EstimateErrorsDialog
from .dialogs import EstimateSNRDialog
from .dialogs import FindContinuumDialog
from .dialogs import FindLinesDialog
from .dialogs import FindSegmentsDialog
from .dialogs import FitContinuumDialog
from .dialogs import FitLinesDialog
from .dialogs import OperateSpectrumDialog
from .dialogs import ResampleSpectrumDialog
from .dialogs import SendSpectrumDialog
from .dialogs import SolverDialog
from .dialogs import SolverEWDialog
from .dialogs import SyntheticSpectrumDialog
from .dialogs import InterpolateSpectrumDialog
from .dialogs import InterpolateSolverDialog
from .dialogs import VelocityProfileDialog

from .CustomizableRegion import CustomizableRegion
from .Spectrum import Spectrum
from .Meter import Meter
from .StatusBar import StatusBar


try:
    from .SAMPManager import SAMPManager
except:
    pass

def resource_path(relative):
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
        # Since we are inside "ispec/gui/", we go up two levels to find "input/"
        basedir = os.path.dirname(basedir[:-1])
        basedir = os.path.dirname(basedir[:-1])
    return os.path.join(basedir, relative)


class iSpecBaseApp(tkinter.Tk):

    def __init_attributes__(self):
        self.velocity_telluric_lower_limit = -100 # km/s
        self.velocity_telluric_upper_limit = 100 # km/s
        self.velocity_telluric_step = 0.5 # km/s
        self.velocity_atomic_lower_limit = -200 # km/s
        self.velocity_atomic_upper_limit = 200 # km/s
        self.velocity_atomic_step = 1.0 # km/s

        self.ccf_mask = {}
        self.ccf_template = {}
        self.atomic_linelist = {}
        self.chemical_elements = None
        self.molecules = None
        self.telluric_linelist = None
        self.solar_abundances = {}

        self.lists = {}
        self.lists['atmospheres'] = self.__get_filelist('input/atmospheres/', 'parameters.tsv')
        self.lists['abundances'] = self.__get_filelist('input/abundances/', 'stdatom.dat')
        self.lists['atomic_lines'] = self.__get_filelist('input/linelists/transitions/', 'atomic_lines.tsv')
        self.lists['masks'] = self.__get_filelist('input/linelists/CCF/', 'mask.lst')
        self.lists['templates'] = self.__get_filelist('input/spectra/templates/', 'template.txt.gz')

        self.default_lists = {}
        self.default_lists['synth_code'] = 0
        self.default_lists['ew_code'] = 0
        self.default_lists['atmospheres'] = 0
        self.default_lists['abundances'] = 0
        self.default_lists['atomic_lines'] = 0
        self.default_lists['grid'] = 0

        # Radiative transfer codes
        self.lists['synth_code'] = []
        self.lists['ew_code'] = []
        if ispec.is_spectrum_support_enabled():
            self.lists['synth_code'].append("SPECTRUM")
            self.lists['ew_code'].append("SPECTRUM")
        if ispec.is_turbospectrum_support_enabled():
            self.lists['synth_code'].append("Turbospectrum")
            self.lists['ew_code'].append("Turbospectrum")
        if ispec.is_sme_support_enabled():
            self.lists['synth_code'].append("SME")
        if ispec.is_moog_support_enabled():
            self.lists['synth_code'].append("MOOG")
            self.lists['ew_code'].append("MOOG")
            self.default_lists['ew_code'] = len(self.lists['ew_code'])-1 # Prefer moog before spectrum
        if ispec.is_width_support_enabled():
            self.lists['ew_code'].append("Width")
            self.default_lists['ew_code'] = len(self.lists['ew_code'])-1 # Prefer width before moog
        if ispec.is_synthe_support_enabled():
            self.lists['synth_code'].append("Synthe")

        self.lists['grid'] = self.__get_filelist('input/grid/', 'parameters.tsv')

        ######
        # Prefered defaults:
        found = np.where(self.lists['atmospheres']['name'] == 'MARCS.GES')[0]
        if len(found) == 1:
            self.default_lists['atmospheres'] = found[0]
        found = np.where(self.lists['abundances']['name'] == 'Grevesse.2007')[0]
        if len(found) == 1:
            self.default_lists['abundances'] = found[0]
        found = np.where(self.lists['atomic_lines']['name'] == 'GESv6_atom_hfs_iso.420_920nm')[0]
        if len(found) == 1:
            self.default_lists['atomic_lines'] = found[0]
        found = np.where(self.lists['grid']['name'] == 'SPECTRUM_MARCS.GES_GESv6_atom_hfs_iso.480_680nm_light')[0]
        if len(found) == 1:
            self.default_lists['grid'] = found[0]
        found = np.where(self.lists['grid']['name'] == 'SPECTRUM_MARCS.GES_GESv6_atom_hfs_iso.480_680nm')[0]
        if len(found) == 1:
            self.default_lists['grid'] = found[0]
        ######

        ######
        # Different way of defining defaults for templates and masks since it is a more difficult case:
        prefered = self.lists['masks']['name'] == 'Narval.Sun.370_1048nm'
        self.lists['masks'] = np.hstack((self.lists['masks'][prefered], self.lists['masks'][~prefered]))
        #
        prefered = self.lists['templates']['name'] == 'Synth.Sun.350_1100nm'
        self.lists['templates'] = np.hstack((self.lists['templates'][prefered], self.lists['templates'][~prefered]))
        ######

        self.velocity_template_lower_limit = -200 # km/s
        self.velocity_template_upper_limit = 200 # km/s
        self.velocity_template_step = 1.0 # km/s
        self.modeled_layers_pack = {} # Synthesize spectrum (atmospheric models)
        self.grid = {} # Interpolate spectrum
        self.find_continuum_regions_wave_step = 0.05
        self.find_continuum_regions_sigma = 0.001
        self.find_continuum_regions_max_continuum_diff = 1.0
        self.find_lines_min_depth = 0.05 # (% of the continuum)
        self.find_lines_max_depth = 1.00 # (% of the continuum)
        self.show_errors = False


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
        # Base reference from where to update the spectra list
        self.menu_active_spectrum_base = 0

    def __periodic_queue_check(self):
        # 100 ms
        while self.queue.qsize():
            try:
                command = self.queue.get(0)
                if len(command) == 3:
                    # func, args, kargs tuple
                    command[0](*command[1], **command[2])
                else:
                    # string command (usefull when called from other processes
                    # such as when a synthetic spectrum is generated)
                    eval(command)
                self.queue.task_done()
            except Empty:
                pass
        self.after(100, self.__periodic_queue_check)

    def __init__(self, spectra, regions, filenames):
        tkinter.Tk.__init__(self)
        self.protocol('WM_DELETE_WINDOW', self.on_close)

        # Window icon
        img = tkinter.PhotoImage(file=resource_path("images/iSpec.gif"))
        self.tk.call('wm', 'iconphoto', self._w, img)
        #self.iconbitmap(bitmap="@"+resource_path("images/iSpec.xbm")) # Black and white

        self.queue = Queue()
        # Start the periodic call in the GUI to check if the queue contains
        # anything
        self.__periodic_queue_check()

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


        if regions is None:
            continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])
            lines = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])
            segments = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])

            self.regions = {}
            self.regions['continuum'] = continuum
            self.regions['lines'] = lines
            self.regions['segments'] = segments
        else:
            self.regions = regions

        if filenames is None:
            self.filenames = {}
            self.filenames['continuum'] = None
            self.filenames['lines'] = None
            self.filenames['fitted_lines'] = None
            self.filenames['segments'] = None
        else:
            self.filenames = {}
            self.filenames['continuum'] = filenames['continuum']
            self.filenames['lines'] = filenames['lines']
            self.filenames['fitted_lines'] = None
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


    def __get_filelist(self, dirname, match):
        import os
        import glob
        ispec_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

        filelist = []
        for root in glob.glob(os.path.join(os.path.join(ispec_dir, dirname), "*")):
            if os.path.isdir(root) and os.path.exists(os.path.join(root, match)):
                filelist.append((os.path.basename(root), resource_path(os.path.join(root, match))))
        filelist = np.array(filelist, dtype=[('name', '|U100'), ('path', '|U500')])
        filelist.sort(order=['name'])
        return filelist


    def create_main_window(self):
        self.create_window()
        self.create_menu()
        self.create_plotting_area()
        self.create_controls()
        self.create_statusbar()

    def create_window(self ):
        self.frame = tkinter.Frame(self)
        self.frame.pack(fill=tkinter.BOTH, expand=1)
        self.wm_title("iSpec")

    def create_menu(self):
        # create a menu
        menu = tkinter.Menu(self)
        self.config(menu=menu)

        filemenu = tkinter.Menu(menu)
        menu.add_cascade(label="Files", menu=filemenu)
        filemenu.add_command(label="Open spectra", command=self.on_open_spectra)
        filemenu.add_command(label="Open continuum regions", command=self.on_open_continuum)
        filemenu.add_command(label="Open line regions", command=self.on_open_lines)
        filemenu.add_command(label="Open segment regions", command=self.on_open_segments)
        filemenu.add_separator()
        filemenu.add_command(label="Save plot image as...", command=self.on_save_plot)
        filemenu.add_command(label="Save spectrum as...", command=self.on_save_spectrum)
        self.spectrum_function_items.append((filemenu, filemenu.entrycget(tkinter.END, "label")))
        filemenu.add_command(label="Save continuum regions as...", command=self.on_save_continuum_regions)
        filemenu.add_command(label="Save line regions as...", command=self.on_save_line_regions)
        filemenu.add_command(label="Save segment as...", command=self.on_save_segments)
        filemenu.add_separator()
        filemenu.add_command(label="Export fitted lines as...", command=self.on_save_fitted_line_regions)
        filemenu.add_separator()
        filemenu.add_command(label="Exit", command=self.on_close)

        operationmenu = tkinter.Menu(menu)
        menu.add_cascade(label="Operations", menu=operationmenu)

        continuummenu = tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Fit continuum...", menu=continuummenu)
        continuummenu.add_command(label="Splines", command=self.on_fit_continuum)
        self.spectrum_function_items.append((continuummenu, continuummenu.entrycget(tkinter.END, "label")))
        continuummenu.add_command(label="Polynomy", command=self.on_fit_continuum_polynomy)
        self.spectrum_function_items.append((continuummenu, continuummenu.entrycget(tkinter.END, "label")))
        continuummenu.add_command(label="Template", command=self.on_fit_continuum_template)
        self.spectrum_function_items.append((continuummenu, continuummenu.entrycget(tkinter.END, "label")))
        continuummenu.add_command(label="Fixed value", command=self.on_fit_continuum_fixed_value)
        self.spectrum_function_items.append((continuummenu, continuummenu.entrycget(tkinter.END, "label")))

        if len(self.lists['atomic_lines']) > 0:
            operationmenu.add_command(label="Fit lines", command=self.on_fit_lines)
            self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_separator()

        clearmenu = tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Clear...", menu=clearmenu)
        clearmenu.add_command(label="Fitted continuum", command=self.on_remove_fitted_continuum)
        clearmenu.add_command(label="Fitted lines", command=self.on_remove_fitted_lines)
        clearmenu.add_command(label="Continuum regions", command=self.on_remove_continuum_regions)
        clearmenu.add_command(label="Line regions", command=self.on_remove_line_masks)
        clearmenu.add_command(label="Segments", command=self.on_remove_segments)

        operationmenu.add_separator()
        operationmenu.add_command(label="Find continuum regions", command=self.on_find_continuum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        if len(self.lists['atomic_lines']) > 0:
            operationmenu.add_command(label="Find line masks", command=self.on_find_lines)
            self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Adjust line masks", command=self.on_adjust_lines)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Create segments around line masks", command=self.on_create_segments_around_lines)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))

        operationmenu.add_separator()


        velocitymenu = tkinter.Menu(operationmenu)
        operationmenu.add_cascade(label="Correct velocity relative to...", menu=velocitymenu)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))

        velocitymenu.add_command(label="Atomic line mask (radial velocity)", command=self.on_correct_velocity_atomic)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))
        velocitymenu.add_command(label="Telluric line mask  (barycentric velocity)", command=self.on_correct_velocity_telluric)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))
        velocitymenu.add_command(label="Template", command=self.on_correct_velocity_template)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))

        operationmenu.add_separator()

        operationmenu.add_command(label="Calculate errors based on SNR", command=self.on_calculate_errors)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Add noise to spectrum fluxes", command=self.on_add_noise)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_separator()
        operationmenu.add_command(label="Degrade resolution", command=self.on_degrade_resolution)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Continuum normalization", command=self.on_continuum_normalization)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Clean fluxes and errors", command=self.on_clean_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Clean telluric regions", command=self.on_clean_tellurics)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Wavelength range reduction", command=self.on_cut_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Apply mathematical expression", command=self.on_operate_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Resample spectrum", command=self.on_resample_spectrum)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))
        operationmenu.add_command(label="Combine all spectra", command=self.on_combine_spectra)
        self.spectrum_function_items.append((operationmenu, operationmenu.entrycget(tkinter.END, "label")))


        parametersmenu = tkinter.Menu(menu)
        menu.add_cascade(label="Parameters", menu=parametersmenu)

        velocitymenu = tkinter.Menu(operationmenu)
        parametersmenu.add_cascade(label="Determine velocity relative to...", menu=velocitymenu)
        self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))

        if len(self.lists['masks']) > 1: # More than one because 1 for tellurics and 1 for stellar lines at least
            velocitymenu.add_command(label="Atomic line mask (radial velocity)", command=self.on_determine_velocity_atomic)
            self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))
            velocitymenu.add_command(label="Telluric line mask  (barycentric velocity)", command=self.on_determine_velocity_telluric)
            self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))
        velocitymenu.add_command(label="Template", command=self.on_determine_velocity_template)
        self.spectrum_function_items.append((velocitymenu, velocitymenu.entrycget(tkinter.END, "label")))
        parametersmenu.add_command(label="Calculate barycentric velocity", command=self.on_determine_barycentric_vel)
        parametersmenu.add_separator()
        parametersmenu.add_command(label="Estimate SNR", command=self.on_estimate_snr)
        self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))
        parametersmenu.add_separator()
        if len(self.lists['grid']) > 0:
            parametersmenu.add_command(label="Determine parameters and abundances with grid", command=self.on_determine_parameters_with_grid)
            self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))
        if (ispec.is_spectrum_support_enabled() \
                #or ispec.is_turbospectrum_support_enabled() \
                #or ispec.is_moog_support_enabled() \
                #or ispec.is_width_support_enabled() \
                ) and \
                len(self.lists['atmospheres']) > 0 and len(self.lists['abundances']) > 0 and len(self.lists['atomic_lines']) > 0:
            parametersmenu.add_command(label="Determine parameters and abundances with synthesis", command=self.on_determine_parameters)
            self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))
            parametersmenu.add_command(label="Determine abundances with equivalent widths", command=self.on_determine_abundances_from_ew)
            self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))
            parametersmenu.add_command(label="Determine parameters with equivalent widths", command=self.on_determine_parameters_from_ew)
            self.spectrum_function_items.append((parametersmenu, parametersmenu.entrycget(tkinter.END, "label")))

        self.menu_active_spectrum_num = tkinter.IntVar()
        self.menu_active_spectrum_num.set('1')

        self.menu_active_spectrum = tkinter.Menu(menu)
        menu.add_cascade(label="Spectra", menu=self.menu_active_spectrum)

        self.menu_active_spectrum.add_command(label="Duplicate spectrum", command=self.on_duplicate_spectrum)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.menu_active_spectrum.add_command(label="Close spectrum", command=self.on_close_spectrum)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.menu_active_spectrum.add_command(label="Close all spectra", command=self.on_close_all_spectra)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.menu_active_spectrum.add_command(label="Close all spectra except the active one", command=self.on_close_all_spectra_except_active)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.menu_active_spectrum.add_separator()

        if len(self.lists['grid']) > 0:
                self.menu_active_spectrum.add_command(label="Interpolate spectrum from grid", command=self.on_interpolate)

        if (ispec.is_spectrum_support_enabled() \
                #or ispec.is_turbospectrum_support_enabled() \
                #or ispec.is_sme_support_enabled() \
                #or ispec.is_moog_support_enabled() \
                #or ispec.is_synthe_support_enabled() \
                ) and \
                len(self.lists['atmospheres']) > 0 and len(self.lists['abundances']) > 0 and len(self.lists['atomic_lines']) > 0:
                self.menu_active_spectrum.add_command(label="Synthesize spectrum", command=self.on_synthesize)

        if self.samp_manager is not None:
            self.menu_active_spectrum.add_command(label="Send spectrum to...", command=self.on_send_spectrum)
            self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.show_errors = tkinter.BooleanVar()
        self.menu_active_spectrum.add_checkbutton(label="Show errors in plot", onvalue=True, offvalue=False, variable=self.show_errors, command=self.on_show_errors)
        self.spectrum_function_items.append((self.menu_active_spectrum, self.menu_active_spectrum.entrycget(tkinter.END, "label")))
        self.menu_active_spectrum.add_separator()

        # Base reference from where to update the spectra list
        self.menu_active_spectrum_base = self.menu_active_spectrum.index("Show errors in plot") + 2

        helpmenu = tkinter.Menu(menu)
        menu.add_cascade(label="Help", menu=helpmenu)

        helpmenu.add_command(label="License", command=self.on_license)
        helpmenu.add_command(label="About...", command=self.on_about)
        self.update_menu_active_spectrum()

    def create_plotting_area(self):
        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.plot_frame = tkinter.Frame(self.frame)
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.fig, master=self.plot_frame)
        self.canvas.draw()

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
        self.canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)
        self.plot_frame.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=1)

    def create_controls(self):
        # controls
        self.control_frame = tkinter.Frame(self.frame)

        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

        label = tkinter.Label(self.control_frame, text="Action:")
        label.pack(side=tkinter.LEFT)
        self.action_entry = tkinter.StringVar()
        self.action_entry.set("Stats")
        self.action_buttons = {}
        for text in ["Stats", "Create", "Modify", "Remove"]:
            self.action_buttons[text] = tkinter.Radiobutton(self.control_frame, text=text, variable=self.action_entry, value=text, command=self.on_action_change)
            self.action_buttons[text].pack(side=tkinter.LEFT)

        label = tkinter.Label(self.control_frame, text=" || Element:")
        label.pack(side=tkinter.LEFT)
        self.elements_entry = tkinter.StringVar()
        self.elements_entry.set("Continuum")
        self.elements_buttons = {}
        for text in ["Continuum", "Lines", "Line marks", "Segments"]:
            if text == "Line marks":
                self.elements_buttons[text] = tkinter.Radiobutton(self.control_frame, text=text, variable=self.elements_entry, value=text, command=self.on_element_change, state=tkinter.DISABLED)
            else:
                self.elements_buttons[text] = tkinter.Radiobutton(self.control_frame, text=text, variable=self.elements_entry, value=text, command=self.on_element_change)
            self.elements_buttons[text].pack(side=tkinter.LEFT)

        #self.button = Tkinter.Button(self.control_frame, text="QUIT", fg="red", command=self.on_close)
        #self.button.pack(side=Tkinter.LEFT)

        self.progress_bar = Meter(self.control_frame)
        self.progress_bar.pack(side=tkinter.LEFT)
        self.control_frame.pack()

        frame = tkinter.Frame(self)
        self.stats_scrollbar = tkinter.Scrollbar(frame, orient=tkinter.VERTICAL)
        self.stats = tkinter.Listbox(frame, height=5, yscrollcommand=self.stats_scrollbar.set, font=('courier',10,'normal'), selectmode=tkinter.EXTENDED)
        self.stats_scrollbar.pack(side=tkinter.RIGHT, fill=tkinter.Y)
        self.stats.pack(fill=tkinter.BOTH, expand=1)
        frame.pack(fill=tkinter.X)

    def create_statusbar(self):
        ## create a statusbar
        self.status = StatusBar(self)
        self.status.pack(side=tkinter.BOTTOM, fill=tkinter.X)
        self.status.set("hi!")


    def draw_figure_for_first_time(self):
        """ Redraws the figure
        """
        self.axes.grid(True, which="both")
        self.axes.set_title("Spectra", fontsize="10")
        self.axes.set_xlabel("wavelength (nm)", fontsize="10")
        self.axes.set_ylabel("flux", fontsize="10")

        for spec in self.spectra:
            if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
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
        if isinstance(name, bytes):
            name = name.decode('utf-8')
        num_repeated = 0
        max_num = 0
        for spec in self.spectra:
            if spec.name.startswith(name):
                try:
                    # Does it has already a suffix?
                    num = int(spec.name.split(self.dupiclated_name_separator)[-1])
                    if num > max_num:
                        max_num = num
                except ValueError:
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

        if free_color_found is None:
            good_color = False
            while not good_color:
                # Random number converted to hexadecimal
                random_num = np.random.randint(0, 16777215)
                random_color = "#%x" % random_num
                # Test that the generated color is correct
                try:
                    matplotlib.colors.colorConverter.to_rgba(random_color)
                    good_color = True
                except ValueError:
                    pass

            free_color_found = random_color

        return free_color_found


    def on_license(self):
        license = """The Integrated Spectroscopic Framework is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

the Integrated Spectroscopic Framework is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along with the Integrated Spectroscopic Framework.  If not, see:

www.gnu.org/licenses/"""
        self.info("iSpec License", license)

    def on_about(self):
        description = """iSpec is a tool for the treatment and analysis of high-resolution and high singal-to-noise stellar spectra (mainly FGK stars) developed by Sergi Blanco-Cuaresma.
"""
        if (ispec.is_spectrum_support_enabled() \
                or ispec.is_turbospectrum_support_enabled() \
                or ispec.is_sme_support_enabled() \
                or ispec.is_moog_support_enabled() \
                or ispec.is_width_support_enabled() \
                or ispec.is_synthe_support_enabled()):
            description += """
iSpec uses the following radiative transfer codes:

1) SPECTRUM - Richard O. Gray - Version 2.77
2) Turbospectrum - Bertrand Plez - v19.1
3) MOOG - Chris Sneden - November 2019
4) SYNTHE/WIDTH9 - R. L. Kurucz / Atmos port - 2015
5) SME - Valenti & Piskunov - 574
"""
        self.info("About iSpec", description)



    def on_action_change(self):
        self.enable_elements()
        self.stats.delete(0, tkinter.END) # Delete all

        self.action = self.action_entry.get()
        if self.action in ["Stats", "Create"]:
            self.elements_buttons['Line marks'].config(state=tkinter.DISABLED)
            if self.elements_entry.get() == "Line marks":
                self.elements_entry.set("continuum")
                self.elements = "continuum"
                self.elements_buttons["Continuum"].select()

    def enable_elements(self):
        for text in ["Continuum", "Lines", "Line marks", "Segments"]:
            self.elements_buttons[text].config(state=tkinter.NORMAL)

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
                if self.samp_manager is not None:
                    self.samp_manager.shutdown()
                self.frame.quit() # stops mainloop
                self.frame.destroy()  # this is necessary on Windows to prevent
                                      # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        else:
            if self.samp_manager is not None:
                self.samp_manager.shutdown()
            self.frame.quit() # stops mainloop
            self.frame.destroy()  # this is necessary on Windows to prevent
                                  # Fatal Python Error: PyEval_RestoreThread: NULL tstate



    def on_close_all_spectra_except_active(self):
        self.__close_all_spectra(except_active=True)

    def on_close_all_spectra(self):
        self.__close_all_spectra(except_active=False)

    def __close_all_spectra(self, except_active=True):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        if len(self.spectra) <= 1 and except_active:
            return

        some_not_saved = False
        for spec in self.spectra:
            if except_active and self.active_spectrum == spec:
                continue
            if spec.not_saved:
                some_not_saved = True
                break

        if except_active:
            msg = "Are you sure you want to close ALL the spectra except the active one?"
        else:
            msg = "Are you sure you want to close ALL the spectra?"
        if some_not_saved:
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            title = "Close all the spectra"
            if not self.question(title, msg):
                return

        for spec in self.spectra:
            if except_active and self.active_spectrum == spec:
                continue
            spec.plot_id.remove()
            # Remove errors if they exists
            if spec is not None and spec.errors_plot_id1 is not None:
                spec.errors_plot_id1.remove()
                spec.errors_plot_id1 = None
            if spec is not None and spec.errors_plot_id2 is not None:
                spec.errors_plot_id2.remove()
                spec.errors_plot_id2 = None
            # Remove fitted continuum if it exists
            if spec is not None and spec.continuum_plot_id is not None:
                spec.continuum_plot_id.remove()
                spec.continuum_plot_id = None
                spec.continuum_model = None
                spec.continuum_data = None
            # Remove fitted lines if they exist
            for region in self.region_widgets["lines"]:
                if spec in region.line_model:
                    if region.line_plot_id[spec] is not None:
                        region.line_plot_id[spec].remove()
                    del region.line_plot_id[spec]
                    del region.line_model[spec]
                    del region.line_extra[spec]
        if except_active:
            self.spectra = [self.active_spectrum]
            self.active_spectrum_history = [self.active_spectrum]
        else:
            self.spectra = []
            self.active_spectrum = None
            self.active_spectrum_history = []

        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.update_scale()
        if except_active:
            self.flash_status_message("All spectra closed except the active one.")
        else:
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
            if self.active_spectrum is not None and self.active_spectrum.path is not None:
                filename = self.active_spectrum.path.split('/')[-1]
                filename_length = len(filename)
                dirname = self.active_spectrum.path[:-filename_length]
            elif self.active_spectrum is not None:
                filename = self.active_spectrum.name
                dirname = os.getcwd()
            else:
                filename = ""
                dirname = os.getcwd()
        else:
            if self.filenames[elements] is not None:
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
                answer = tkinter.filedialog.askopenfilenames(title="Open %s..." % elements, initialdir=dirname, filetypes=ftypes, defaultextension=".txt")
                unique_answer = np.unique(answer)
                if len(unique_answer) > 0 and len(unique_answer[0]) > 0:
                    answer_ok = True
                else:
                    answer_ok = False
            else:
                ftypes = [('All files', '*'), ('Plain text', '*.txt')]
                if sys.platform == "darwin":
                    ftypes = [] # Not working properly in MacOSX
                answer = tkinter.filedialog.askopenfilename(title="Open %s..." % elements, initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
                if isinstance(answer, str) and len(answer) > 0:
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
                            if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
                                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
                            new_spectrum_data = ispec.read_spectrum(path)
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
                            self.regions[elements] = ispec.read_continuum_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        elif elements == "lines":
                            self.remove_fitted_lines()
                            self.regions[elements] = ispec.read_line_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        else:
                            # 'segments'
                            self.regions[elements] = ispec.read_segment_regions(path)
                            self.draw_regions(elements)
                            self.not_saved[elements] = False
                            self.update_title()
                        self.flash_status_message("Opened file %s" % path)
                    self.filenames[elements] = path
                    self.canvas.draw()
                    action_ended = True
                except Exception as e:
                    msg = 'A file does not have a compatible format.'
                    title = 'File format incompatible'
                    self.error(title, msg)
                    continue
            else:
                self.flash_status_message("Discarded.")
                action_ended = True



    def on_save_plot(self):
        if self.check_operation_in_progress():
            return
        #file_choices = "PNG (*.png)|*.png"

        if self.active_spectrum is not None:
            if self.active_spectrum.path is not None:
                filename = self.active_spectrum.path.split('/')[-1] + ".png"
                filename_length = len(filename)
                dirname = self.active_spectrum.path[:-filename_length]
            else:
                filename = self.active_spectrum.name.split(self.dupiclated_name_separator)[0] + ".png"
                dirname = os.getcwd()
        else:
            filename = "iSpec_plot_image.png"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            ftypes = [('PNG (*.png)', '*.png')]
            if sys.platform == "darwin":
                ftypes = [] # Not working properly in MacOSX
            answer = tkinter.filedialog.asksaveasfilename(title="Save plot as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".png")
            if isinstance(answer, str) and len(answer) > 0:
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

        if self.active_spectrum.path is not None:
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
            answer = tkinter.filedialog.asksaveasfilename(title="Save spectrum as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
            if isinstance(answer, str) and len(answer) > 0:
                answer_ok = True
            else:
                answer_ok = False

            if answer_ok:
                path = answer
                self.status_message("Saving %s..." % path)
                # Save, compress if the filename ends with ".gz"
                ispec.write_spectrum(self.active_spectrum.data, path)
                self.active_spectrum.not_saved = False
                self.update_title()

                # Change name and path
                self.active_spectrum.path = path
                name = self.get_name(path.split('/')[-1]) # If it already exists, add a suffix
                self.active_spectrum.name = name
                self.active_spectrum.plot_id.set_label("[A] " + self.active_spectrum.name)
                self.update_legend()
                self.canvas.draw()

                # Menu active spectrum
                index = int(self.menu_active_spectrum_num.get())-1
                self.menu_active_spectrum.entryconfig(self.menu_active_spectrum_base+index, label=self.active_spectrum.name)

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
        #file_choices = "All|*"
        saved = False
        elements = elements.lower()

        if elements == "fitted_lines":
            if self.active_spectrum.linemasks is None:
                msg = "There is no fitted lines to be exported"
                title = 'No fitted lines'
                self.error(title, msg)
                return
        else:
            if len(self.region_widgets[elements]) == 0:
                msg = "There is no regions to be saved"
                title = 'Empty regions'
                self.error(title, msg)
                return

        if self.filenames[elements] is not None:
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
            answer = tkinter.filedialog.asksaveasfilename(title="Save regions as...", initialdir=dirname, initialfile=filename, filetypes=ftypes, defaultextension=".txt")
            if isinstance(answer, str) and len(answer) > 0:
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

                if elements == "fitted_lines":
                    ispec.write_line_regions(self.active_spectrum.linemasks, path, extended=True)
                else:
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
                                if region.note is not None:
                                    note = region.note.get_text()
                                else:
                                    note = ""
                                output.write("%.5f" % region.get_wave_peak() + "\t" + "%.5f" % region.get_wave_base() + "\t" + "%.5f" % region.get_wave_top() + "\t" + note + "\n")
                            else:
                                output.write("%.5f" % region.get_wave_base() + "\t" + "%.5f" % region.get_wave_top() + "\n")

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

    def on_save_fitted_line_regions(self):
        if self.check_operation_in_progress():
            return
        self.save_regions("fitted_lines")

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
        self.active_spectrum.plot_id.remove()

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
        if self.active_spectrum is not None and self.active_spectrum.errors_plot_id1 is not None:
            self.active_spectrum.errors_plot_id1.remove()
            self.active_spectrum.errors_plot_id1 = None
        if self.active_spectrum is not None and self.active_spectrum.errors_plot_id2 is not None:
            self.active_spectrum.errors_plot_id2.remove()
            self.active_spectrum.errors_plot_id2 = None


    def on_motion(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes is None: return
        if event.inaxes.get_navigate_mode() is not None: return
        if self.programmed_flash_status is not None: return
        if self.operation_in_progress: return
        self.status_message("Cursor on wavelength %.4f" % event.xdata + " and flux %.4f" % event.ydata)

    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes is None: return
        if event.inaxes.get_navigate_mode() is not None: return
        if self.check_operation_in_progress():
            return

        new_halfwidth = 0.05
        # If the left button is clicked when action "create" is active,
        # create a new region of the selected element type
        # NOTE: line marks only can be created from a line region
        if self.action == "Create" and self.subelements != "marks" and event.button == 1 and event.key is None:
            if self.elements == "continuum":
                region = CustomizableRegion(self, "continuum", event.xdata - new_halfwidth, event.xdata + new_halfwidth)
                region.connect()
                self.region_widgets['continuum'].append(region)
            elif self.elements == "lines":
                note_text = self.ask_value('Note for the new line region:', 'Note', '')
                if note_text is None:
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
        self.menu_active_spectrum.delete(self.menu_active_spectrum_base, tkinter.END)

        if len(self.spectra) == 0:
            # No spectra loaded
            self.menu_active_spectrum.add_radiobutton(label="None", variable=self.menu_active_spectrum_num, value=str(1), indicatoron=0, state=tkinter.DISABLED)
        else:
            # Add as many options as spectra
            for i in np.arange(len(self.spectra)):
                if self.spectra[i] is None:
                    continue

                self.menu_active_spectrum.add_radiobutton(label=self.spectra[i].name, variable=self.menu_active_spectrum_num, value=str(i+1), indicatoron=1, command=self.on_change_active_spectrum)

                if self.active_spectrum == self.spectra[i]:
                    self.menu_active_spectrum_num.set(str(i+1))

        if len(self.spectra) > 0:
            for menu, label in self.spectrum_function_items:
                index = menu.index(label)
                menu.entryconfig(index, state=tkinter.NORMAL)
        else:
            for menu, label in self.spectrum_function_items:
                index = menu.index(label)
                menu.entryconfig(index, state=tkinter.DISABLED)

    def draw_active_spectrum(self):
        if self.active_spectrum is not None:
            # Remove spectrum plot if exists
            if self.active_spectrum.plot_id is not None:
                self.active_spectrum.plot_id.remove()

            # zorder = 1, always in the background
            self.active_spectrum.plot_id = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'], lw=1, color=self.active_spectrum.color, linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor=self.active_spectrum.color, zorder=1, label="[A] "+self.active_spectrum.name)[0]

            # Draw errors
            if self.show_errors.get():
                if self.active_spectrum.errors_plot_id1 is not None:
                    self.active_spectrum.errors_plot_id1.remove()
                if self.active_spectrum.errors_plot_id2 is not None:
                    self.active_spectrum.errors_plot_id2.remove()

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
        self.stats.insert(tkinter.END, "%-40s: %s" % (str(k), str(v)))

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
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
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
        if self.check_operation_in_progress():
            return
        self.queue.put((self.flash_status_message, ["Receiving spectrum from external application..."], {'progress':False}))
        self.operation_in_progress = True
        self.queue.put((self.on_receive_spectrum_thread, [new_spectrum_data, name], {}))

    def on_receive_spectrum_thread(self, new_spectrum_data, name):
        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
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
        if self.samp_manager is None:
            return

        if not self.samp_manager.is_connected():
            msg = "No compatible external application can be detected because iSpec is not connected to any SAMP hub.\n\n* A SAMP hub can be created by using TOPCAT application."
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
        if key not in self.dialog:
            self.dialog[key] = SendSpectrumDialog(self, "Send spectrum to external application", names)
        self.dialog[key].show()

        if self.dialog[key].results is None:
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

        if selected_application_index is None:
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
        self.queue.put((self.on_send_spectrum_finnish, [application, errors], {}))

    def on_send_spectrum_finnish(self, application, errors):
        if errors:
            self.flash_status_message("Communication failed with %s." % (application))
        else:
            self.flash_status_message("Spectrum sent to %s." % (application))
        self.operation_in_progress = False


    def remove_drawn_errors_spectra(self):
        for spec in self.spectra:
            if spec.errors_plot_id1 is not None:
                spec.errors_plot_id1.remove()
                spec.errors_plot_id1 = None
            if spec.errors_plot_id2 is not None:
                spec.errors_plot_id2.remove()
                spec.errors_plot_id2 = None
        self.canvas.draw()

    def draw_errors_spectra(self):
        for spec in self.spectra:
            # Remove continuum plot if exists
            if spec.errors_plot_id1 is not None:
                spec.errors_plot_id1.remove()
            if spec.errors_plot_id2 is not None:
                spec.errors_plot_id2.remove()

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
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|U100')])
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
        if self.active_spectrum is None or len(self.active_spectrum.data['waveobs']) == 0:
            msg = "There is no active spectrum"
            title = "Spectrum not present"
            self.error(title, msg)
            return False
        return True

    def check_continuum_model_exists(self):
        if self.active_spectrum.continuum_model is None:
            msg = "Please, execute a general continuum fit first"
            title = "Continuum model not fitted"
            self.error(title, msg)
            return False
        return True

    def error(self, title, msg):
        tkinter.messagebox.showerror(title, msg)

    def info(self, title, msg):
        tkinter.messagebox.showinfo(title, msg)

    def question(self, title, msg):
        answer_yes = tkinter.messagebox.askquestion(title, msg) == 'yes'
        if not answer_yes:
            self.flash_status_message("Discarded")
        return answer_yes

    def ask_value(self, text, title, default_value):
        response = tkinter.simpledialog.askstring(title, text, initialvalue=default_value)
        return response

    def update_progress(self, value):
        normalized_value = np.min([1, np.max([0, 1.*value / 100])])
        try:
            self.progress_bar.set(normalized_value)
        except:
            # Sometimes this set fails because unknown reasons, we can just ignore the exception
            pass
            #import ipdb
            #ipdb.set_trace()

    def update_title(self):
        title = "iSpec"

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

        self.stats.delete(0, tkinter.END) # Delete all

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

        if region.element_type == "lines" and self.active_spectrum in region.line_plot_id and region.line_model[self.active_spectrum] is not None:
            self.add_stats("Gaussian mean (mu)", "%.4f" % region.line_model[self.active_spectrum].mu())
            self.add_stats("Gaussian mean (mu) error", "%.4f" % region.line_model[self.active_spectrum].emu())
            self.add_stats("Gaussian amplitude (A)", "%.4f" % region.line_model[self.active_spectrum].A())
            self.add_stats("Gaussian standard deviation (sigma)", "%.4f" % region.line_model[self.active_spectrum].sig())
            self.add_stats("Gaussian base level (mean continuum)", "%.4f" % region.line_model[self.active_spectrum].baseline())
            rms = region.line_model[self.active_spectrum].rms
            self.add_stats("Gaussian fit root mean squeare (RMS)", "%.4f" % rms)

            if type(region.line_model[self.active_spectrum]) is ispec.GaussianModel:
                # If it is a gaussian we can directly use a formule (but not if it is a voigt!)
                A = region.line_model[self.active_spectrum].A()
                sig = region.line_model[self.active_spectrum].sig()
                ew = -1.*A*np.sqrt(2*np.pi*sig**2) # nm
                integrated_flux = ew/ region.line_model[self.active_spectrum].baseline() # nm^2
            else:
                integrated_flux = -1 * region.line_model[self.active_spectrum].integrate()
                ew = integrated_flux/ region.line_model[self.active_spectrum].baseline()
            ew = ew * 10000 # From nm to mA
            self.add_stats("Gaussian fit Equivalent Width (EW)", "%.2f" % ew)

        if region.element_type == "lines" and self.active_spectrum in region.line_extra and region.line_extra[self.active_spectrum] is not None:
            # Extras (all in string format separated by ;)
            atomic_wave_peak, species, lower_state, upper_state, loggf, fudge_factor, transition_type, rad, stark, waals, ew, element, telluric_wave_peak, telluric_depth = region.line_extra[self.active_spectrum].split(";")
            self.add_stats("Atomic element", element)
            self.add_stats("Atomic line wavelength", atomic_wave_peak)
            self.add_stats("Atomic lower state (cm^-1)", lower_state)
            self.add_stats("Atomic upper state (cm^-1)", upper_state)
            self.add_stats("Atomic log(gf)", loggf)
            self.add_stats("Atomic radiative damping constant", rad)
            self.add_stats("Atomic stark damping constant", stark)
            self.add_stats("Atomic van der Waals damping constant", waals)
            if telluric_wave_peak != "" and float(telluric_wave_peak) != 0:
                self.add_stats("Tellurics: possibly affected by line at (nm)", telluric_wave_peak)
                self.add_stats("Tellurics: typical line depth", "%.4f" % float(telluric_depth))
            else:
                self.add_stats("Tellurics:", "Unaffected")

        if self.active_spectrum.continuum_model is not None:
            if num_points > 0:
                mean_continuum = np.mean(self.active_spectrum.continuum_model(spectrum_window['waveobs']))
                self.add_stats("Continuum mean for the region", "%.4f" % mean_continuum)
            try:
                residuals = np.abs(self.active_spectrum.continuum_model.residuals())
                rms = np.sqrt(np.sum(np.power(residuals,2))/len(residuals))
                self.add_stats("Continuum fit root mean square (RMS)", "%.4f" % rms)
            except AttributeError:
                # Only spline continuum models have "residuals" method
                pass



    def flash_status_message(self, msg, flash_len_ms=3000, progress=True):
        if self.programmed_flash_status is not None:
            self.after_cancel(self.programmed_flash_status)
            self.programmed_flash_status = None

        self.status.set(msg)
        if progress:
            self.update_progress(100)

        self.programmed_flash_status = self.after(flash_len_ms, self.on_flash_status_off)

    def status_message(self, msg):
        if self.programmed_flash_status is not None:
            self.after_cancel(self.programmed_flash_status)
            self.programmed_flash_status = None
            # Progress bar to zero
            self.update_progress(0)

        self.status.set(msg)

    def remove_drawn_continuum_spectrum(self):
        if self.active_spectrum is not None and self.active_spectrum.continuum_plot_id is not None:
            self.active_spectrum.continuum_plot_id.remove()
            self.active_spectrum.continuum_plot_id = None
            #self.active_spectrum.continuum_model = None
            #self.active_spectrum.continuum_data = None
            self.canvas.draw()

    def remove_drawn_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if self.active_spectrum in region.line_plot_id and region.line_plot_id[self.active_spectrum] is not None:
                region.line_plot_id[self.active_spectrum].remove()
                region.line_plot_id[self.active_spectrum] = None
        self.canvas.draw()

    def remove_fitted_lines(self):
        self.remove_drawn_fitted_lines()
        # Remove fitted lines if they exist
        for region in self.region_widgets["lines"]:
            if self.active_spectrum in region.line_model:
                del region.line_plot_id[self.active_spectrum]
                del region.line_model[self.active_spectrum]
                del region.line_extra[self.active_spectrum]
        if self.active_spectrum is not None:
            self.active_spectrum.linemasks = None
            self.active_spectrum.abundances = None

    def draw_continuum_spectrum(self):
        if self.active_spectrum is None:
            return

        # Remove continuum plot if exists
        if self.active_spectrum.continuum_plot_id is not None:
            self.active_spectrum.continuum_plot_id.remove()

        # zorder = 1, always in the background
        if self.active_spectrum.continuum_data is not None:
            self.active_spectrum.continuum_plot_id = self.axes.plot(self.active_spectrum.continuum_data['waveobs'], self.active_spectrum.continuum_data['flux'], lw=1, color='green', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]

    def draw_fitted_lines(self):
        for region in self.region_widgets["lines"]:
            if self.active_spectrum in region.line_model and region.line_model[self.active_spectrum] is not None:
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


    def on_fit_continuum_polynomy(self):
        self.on_fit_continuum(fit_type="Polynomy")
    def on_fit_continuum_template(self):
        self.on_fit_continuum(fit_type="Template")
    def on_fit_continuum_fixed_value(self):
        self.on_fit_continuum(fit_type="Fixed value")

    def on_fit_continuum(self, fit_type="Splines"):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        # Give priority to the resolution determined by the velocity profile relative to telluric lines
        if self.active_spectrum.resolution_telluric == 0.0:
            R = self.active_spectrum.resolution_atomic
        else:
            R = self.active_spectrum.resolution_telluric

        templates = self.lists['templates']['name'].tolist()
        templates = ["i:"+t for t in templates]
        # Add as many options as spectra
        for i in np.arange(len(self.spectra)):
            if self.spectra[i] is None:
                continue
            templates.append(self.spectra[i].name)

        key = "FitContinuumDialog"+fit_type
        # Initial recommendation: 1 knot every 5 nm
        nknots = np.max([1, int((np.max(self.active_spectrum.data['waveobs']) - np.min(self.active_spectrum.data['waveobs'])) / 5.)])
        degree = 2
        if key not in self.active_spectrum.dialog:
            if fit_type == "Template":
                median_wave_range=5.0
            else:
                median_wave_range=0.05
            max_wave_range=1.0
            #median_wave_range=0.5
            #max_wave_range=0.1
            strong_line_probability = 0.50
            self.active_spectrum.dialog[key] = FitContinuumDialog(self, "Properties for fitting continuum", R, nknots, degree, median_wave_range, max_wave_range, strong_line_probability, templates, fit_type)
        self.active_spectrum.dialog[key].show(suggested_nknots=nknots, updated_templates=templates)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        model = self.active_spectrum.dialog[key].results["Fitting model"]
        if model != "Fixed value":
            fixed_value = None
            if model == "Splines":
                nknots = self.active_spectrum.dialog[key].results["Number of splines"]
            else:
                nknots = None # In Template mode, used to apply a gaussian filter and the default estimation is good enough
            median_wave_range = self.active_spectrum.dialog[key].results["Wavelength step for median selection"]
            if model in ['Splines', 'Polynomy']:
                degree = self.active_spectrum.dialog[key].results["Degree"]
                max_wave_range = self.active_spectrum.dialog[key].results["Wavelength step for max selection"]
                order = self.active_spectrum.dialog[key].results["Filtering order"]
                use_errors_for_fitting = self.active_spectrum.dialog[key].results["Use spectrum's errors as weights for the fitting process"] == 1
                automatic_strong_line_detection = self.active_spectrum.dialog[key].results["Automatically find and ignore strong lines"] == 1
                strong_line_probability = self.active_spectrum.dialog[key].results["Strong line probability threshold"]
            else:
                degree = 0
                max_wave_range = median_wave_range + 0.001
                order='median+max'
                use_errors_for_fitting = True
                automatic_strong_line_detection = True
                strong_line_probability = 0

            R = self.active_spectrum.dialog[key].results["Resolution"]
            if R <= 0:
                R = None
            else:
                self.active_spectrum.resolution_telluric = R
                self.active_spectrum.resolution_atomic = R
            in_continuum = self.active_spectrum.dialog[key].results["Consider only continuum regions"] == 1
            ignore_lines = self.active_spectrum.dialog[key].results["Ignore line regions"] == 1
            each_segment = self.active_spectrum.dialog[key].results["Treat each segment independently"] == 1
            if model == "Template":
                template = self.active_spectrum.dialog[key].results["Use as a template"]
            else:
                template = None

            if median_wave_range < 0 or max_wave_range < 0:
                self.flash_status_message("Bad value.")
                return

            if order == "max+median" and median_wave_range <= max_wave_range:
                msg = "For 'max+median' order, median_wave_range should be greater than max_wave_rage!"
                title = "Bad parameters"
                self.error(title, msg)
                self.flash_status_message("Bad parameters.")
                return
            if order == "median+max" and median_wave_range >= max_wave_range:
                msg = "For 'max+median' order, median_wave_range should be smaller than max_wave_rage!"
                title = "Bad parameters"
                self.error(title, msg)
                self.flash_status_message("Bad parameters.")
                return

            if in_continuum and len(self.region_widgets["continuum"]) == 0:
                msg = 'There are no continuum regions.'
                title = "No continuum regions"
                self.error(title, msg)
                self.flash_status_message("No continuum regions found.")
                return
            if ignore_lines and len(self.region_widgets["lines"]) == 0:
                msg = 'There are no line regions.'
                title = "No line regions"
                self.error(title, msg)
                self.flash_status_message("No line regions found.")
                return
            if each_segment and len(self.region_widgets["segments"]) == 0:
                msg = 'There are no segments.'
                title = "No segments"
                self.error(title, msg)
                self.flash_status_message("No segments found.")
                return
        else:
            fixed_value = self.active_spectrum.dialog[key].results["Fixed value"]
            nknots = None
            median_wave_range = None
            max_wave_range = None
            in_continuum = False
            ignore_lines = False
            each_segment = False
            use_errors_for_fitting = False
            order='median+max'
            automatic_strong_line_detection = False
            strong_line_probability = 0
            template = None
        self.active_spectrum.dialog[key].destroy()


        self.operation_in_progress = True
        self.status_message("Fitting continuum...")
        self.update_progress(10)
        thread = threading.Thread(target=self.on_fit_continuum_thread, args=(nknots,), kwargs={'ignore_lines':ignore_lines, 'in_continuum':in_continuum, 'each_segment': each_segment, 'median_wave_range':median_wave_range, 'max_wave_range':max_wave_range, 'fixed_value':fixed_value, 'model':model, 'degree':degree, 'R':R, 'order':order, 'use_errors_for_fitting': use_errors_for_fitting, 'automatic_strong_line_detection':automatic_strong_line_detection, 'strong_line_probability': strong_line_probability, 'template': template})
        thread.setDaemon(True)
        thread.start()

    def on_fit_continuum_thread(self, nknots, ignore_lines=False, in_continuum=False, each_segment=False, median_wave_range=0.1, max_wave_range=1, fixed_value=None, model="Polynomy", degree=3, R=None, order='median+max', use_errors_for_fitting=True, automatic_strong_line_detection=True, strong_line_probability=0.50, template=None):
        try:
            if each_segment:
                self.__update_numpy_arrays_from_widgets("segments")
                independent_regions = self.regions["segments"]
            else:
                independent_regions = None

            if in_continuum:
                self.__update_numpy_arrays_from_widgets("continuum")
                continuum_regions = self.regions["continuum"]
            else:
                continuum_regions = None

            if ignore_lines:
                self.__update_numpy_arrays_from_widgets("lines")
                ignore_lines = self.regions["lines"]
            else:
                ignore_lines = None

            if template is not None:
                if template.startswith("i:"):
                    # Internal template (solar type)
                    if not template in list(self.ccf_template.keys()):
                        i = np.where(self.lists['templates']['name'] == template[2:])[0][0]
                        self.ccf_template[template] = ispec.read_spectrum(self.lists['templates']['path'][i])
                    template_spectrum = self.ccf_template[template]
                else:
                    # Search template to be used by its name
                    for i in np.arange(len(self.spectra)):
                        if self.spectra[i] is None:
                            continue
                        if self.spectra[i].name == template:
                            template_spectrum = self.spectra[i].data
                            break
            else:
                template_spectrum = None

            self.active_spectrum.continuum_model = ispec.fit_continuum(self.active_spectrum.data, from_resolution=R, independent_regions=independent_regions, continuum_regions=continuum_regions, ignore=ignore_lines, nknots=nknots, degree=degree, median_wave_range=median_wave_range, max_wave_range=max_wave_range, fixed_value=fixed_value, model=model, order=order, automatic_strong_line_detection=automatic_strong_line_detection, strong_line_probability=strong_line_probability, use_errors_for_fitting=use_errors_for_fitting, template=template_spectrum)
            waveobs = self.active_spectrum.data['waveobs']
            self.active_spectrum.continuum_data = ispec.create_spectrum_structure(waveobs, self.active_spectrum.continuum_model(waveobs))

            self.queue.put((self.on_fit_continuum_finish, [nknots], {}))
        except Exception:
            self.operation_in_progress = False
            self.queue.put((self.flash_status_message, ["Not possible to fit continuum with the selected parameters."], {}))
            raise


    def on_fit_continuum_finish(self, nknots):
        self.draw_continuum_spectrum()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Continuum fitted.")

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
            velocity_lower_limit = self.velocity_atomic_lower_limit
            velocity_upper_limit = self.velocity_atomic_upper_limit
            velocity_step = self.velocity_atomic_step
            templates = []
            masks = self.lists['masks']['name'][self.lists['masks']['name'] != 'Synth.Tellurics.500_1100nm']
            mask_size = 2.0
            mask_depth = 0.1
        elif relative_to_telluric_data:
            title = "Velocity profile relative to telluric lines"
            velocity_lower_limit = self.velocity_telluric_lower_limit
            velocity_upper_limit = self.velocity_telluric_upper_limit
            velocity_step = self.velocity_telluric_step
            templates = []
            masks = self.lists['masks']['name'][self.lists['masks']['name'] == 'Synth.Tellurics.500_1100nm']
            mask_size = 2.0
            mask_depth = 0.01
        elif relative_to_template:
            title = "Velocity profile relative to template"
            velocity_lower_limit = self.velocity_template_lower_limit
            velocity_upper_limit = self.velocity_template_upper_limit
            velocity_step = self.velocity_template_step
            masks = []
            mask_size = None
            mask_depth = None
            templates = self.lists['templates']['name'].tolist()
            templates = ["i:"+t for t in templates]
            # Add as many options as spectra
            for i in np.arange(len(self.spectra)):
                if self.spectra[i] is None:
                    continue
                templates.append(self.spectra[i].name)
        else:
            raise Exception("Velocity should be determined relative to something!")

        key = "VelocityProfileDialog:"+str(relative_to_atomic_data)+str(relative_to_telluric_data)+str(relative_to_template)
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = VelocityProfileDialog(self, title, rv_lower_limit=velocity_lower_limit, rv_upper_limit=velocity_upper_limit, rv_step=velocity_step, templates=templates, masks=masks, mask_size=mask_size, mask_depth=mask_depth)
            self.active_spectrum.dialog[key].show()
        elif show_previous_results:
            if relative_to_template:
                self.active_spectrum.dialog[key].show(updated_templates=templates)
            else:
                self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        rv_lower_limit = self.active_spectrum.dialog[key].results["Velocity lower limit (km/s)"]
        rv_upper_limit = self.active_spectrum.dialog[key].results["Velocity upper limit (km/s)"]
        rv_step = self.active_spectrum.dialog[key].results["Velocity steps (km/s)"]
        fourier = self.active_spectrum.dialog[key].results["CCF in Fourier space"] == 1
        model = self.active_spectrum.dialog[key].results["Fitting model"]
        peak_probability = self.active_spectrum.dialog[key].results["Peak probability"]
        if relative_to_template:
            template = self.active_spectrum.dialog[key].results["Cross-correlate with"]
            mask_name = None
            mask_size = None
            mask_depth = None
        else:
            template = None
            if relative_to_atomic_data:
                mask_name = self.active_spectrum.dialog[key].results["Mask linelist"]
            else:
                mask_name = masks[0]
            mask_size = self.active_spectrum.dialog[key].results["Mask size (km/s)"]
            mask_depth = self.active_spectrum.dialog[key].results["Minimum depth"]
        self.active_spectrum.dialog[key].destroy()

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
        thread = threading.Thread(target=self.on_determine_velocity_thread, args=(relative_to_atomic_data, relative_to_telluric_data, relative_to_template, rv_lower_limit, rv_upper_limit, rv_step, template, mask_name, mask_size, mask_depth, fourier, model, peak_probability))
        thread.setDaemon(True)
        thread.start()

    def __filter_telluric_lines(self, telluric_linelist, spectrum, velocity_lower_limit, velocity_upper_limit):
        """
        Select telluric lines inside a given velocity limits (km/s).
        """
        # Light speed in vacuum
        c = 299792458.0 # m/s

        ## Select telluric lines of interest
        # Limit to region of interest
        wmin = spectrum['waveobs'][0]
        wmax = spectrum['waveobs'][-1]
        delta_wmin = wmin * (velocity_lower_limit/ (c/1000.0))
        delta_wmax = wmax * (velocity_upper_limit/ (c/1000.0))
        wfilter = (telluric_linelist['wave_peak'] <= wmax + delta_wmax) & (telluric_linelist['wave_peak'] >= wmin + delta_wmin)
        linelist = telluric_linelist[wfilter]
        # Discard not fitted lines
        rfilter = linelist['rms'] == 9999
        linelist = linelist[~rfilter]
        # Discard too deep or too small lines
        rfilter = (linelist['depth'] <= 0.9) & (linelist['depth'] >= 0.01)
        linelist = linelist[rfilter]
        # Discard outliers FWHM in km/s (which is not wavelength dependent)
        telluric_fwhm = (c/ (linelist['wave_peak']/ linelist['fwhm'])) / 1000.0 # km/s
        fwhm_selected, fwhm_selected_filter = ispec.sigma_clipping(telluric_fwhm, meanfunc=np.median)
        linelist = linelist[fwhm_selected_filter]
        return linelist


    def on_determine_velocity_thread(self, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, velocity_lower_limit, velocity_upper_limit, velocity_step, template, mask_name, mask_size, mask_depth, fourier, model, peak_probability):
        self.queue.put((self.status_message, ["Determining velocity..."], {}))

        if relative_to_atomic_data or relative_to_telluric_data:
            if not mask_name in list(self.ccf_mask.keys()):
                i = np.where(self.lists['masks']['name'] == mask_name)[0][0]
                self.ccf_mask[mask_name] = ispec.read_cross_correlation_mask(self.lists['masks']['path'][i])
            mask_linelist = self.ccf_mask[mask_name]
        elif relative_to_template:
            mask_linelist = None
        else:
            raise Exception("Velocity should be determined relative to something!")

        if relative_to_template:
            if template.startswith("i:"):
                # Internal template (solar type)
                if not template in list(self.ccf_template.keys()):
                    i = np.where(self.lists['templates']['name'] == template[2:])[0][0]
                    self.ccf_template[template] = ispec.read_spectrum(self.lists['templates']['path'][i])
                template_spectrum = self.ccf_template[template]
            else:
                # Search template to be used by its name
                for i in np.arange(len(self.spectra)):
                    if self.spectra[i] is None:
                        continue
                    if self.spectra[i].name == template:
                        template_spectrum = self.spectra[i].data
                        break
            models, ccf = ispec.cross_correlate_with_template(self.active_spectrum.data, template_spectrum, \
                                    lower_velocity_limit=velocity_lower_limit, upper_velocity_limit=velocity_upper_limit, \
                                    velocity_step=velocity_step, \
                                    fourier = fourier, \
                                    model = model, \
                                    only_one_peak = relative_to_telluric_data,
                                    peak_probability = peak_probability)

        else:
            models, ccf = ispec.cross_correlate_with_mask(self.active_spectrum.data, mask_linelist, \
                                    lower_velocity_limit=velocity_lower_limit, upper_velocity_limit=velocity_upper_limit, \
                                    velocity_step=velocity_step, mask_size=mask_size, mask_depth=mask_depth, \
                                    fourier = fourier, \
                                    model = model, \
                                    only_one_peak = relative_to_telluric_data,
                                    peak_probability = peak_probability)

        self.queue.put((self.on_determine_velocity_finish, [models, ccf, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, mask_linelist, model], {}))

    def on_determine_velocity_finish(self, models, ccf, relative_to_atomic_data, relative_to_telluric_data, relative_to_template, mask_linelist, model):

        if len(models) == 0:
            fwhm = 0.0
            telluric_fwhm = 0.0
            R = 0.0
            velocity = 0.0
            self.flash_status_message("Velocity could not be determined!")
        else:
            # Resolution
            c = 299792458.0 # m/s
            fwhm = models[0].fwhm()[0] # km/s (because xcoord is already velocity)
            if relative_to_atomic_data or relative_to_template:
                telluric_fwhm = 0.0
                R = np.int(c/(1000.0*fwhm))
            else:
                # If telluric lines have been used, we can substract its natural FWHM
                # so that we get the real resolution of the instrument (based on the difference in FWHM)
                c = 299792458.0 # m/s
                telluric_fwhm = np.mean((c/ (mask_linelist['wave_peak']/ mask_linelist['fwhm'])) / 1000.0) # km/s
                diff = np.round(fwhm - telluric_fwhm, 2)
                if diff > 0:
                    R = np.int(c/(1000.0*diff))
                else:
                    R = 0
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
        self.active_spectrum.dialog[key].register(ccf, models, telluric_fwhm=telluric_fwhm)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
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

        if self.active_spectrum.resolution_telluric == 0.0:
            resolution = self.active_spectrum.resolution_atomic
        else:
            resolution = self.active_spectrum.resolution_telluric

        key = "FitLinesDialog"
        vel_telluric = self.active_spectrum.velocity_telluric
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = FitLinesDialog(self, "Fit lines", resolution, vel_telluric, self.lists, self.default_lists)
        self.active_spectrum.dialog[key].show(updated_vel_telluric=vel_telluric)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        resolution = self.active_spectrum.dialog[key].results["Resolution"]
        vel_telluric = self.active_spectrum.dialog[key].results["Velocity respect to telluric lines (km/s)"]
        selected_linelist = self.active_spectrum.dialog[key].results["Line list"]
        ### Find filename
        i = np.where(self.lists['atomic_lines']['name'] == selected_linelist)
        atomic_linelist_file = self.lists['atomic_lines']['path'][i][0]
        ####
        free_mu = self.active_spectrum.dialog[key].results["Allow peak position adjustment"] == 1
        check_derivatives = self.active_spectrum.dialog[key].results["Check derivatives before fitting"] == 1
        max_atomic_wave_diff = self.active_spectrum.dialog[key].results["Maximum atomic wavelength difference"]
        self.active_spectrum.dialog[key].destroy()

        telluric_linelist_file = resource_path("input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")

        # Read
        if self.molecules is None:
            self.molecules = ispec.read_molecular_symbols(molecules_file)
        if self.chemical_elements is None:
            self.chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
        if not selected_linelist in list(self.atomic_linelist.keys()):
            self.atomic_linelist[selected_linelist] = ispec.read_atomic_linelist(atomic_linelist_file)
            logging.warning("Limiting linelist to lines that have at least 0.01 depth in the Sun")
            solar = self.atomic_linelist[selected_linelist]['theoretical_depth'] >= 0.01
            self.atomic_linelist[selected_linelist] = self.atomic_linelist[selected_linelist][solar]

        if self.telluric_linelist is None:
            self.telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

        if vel_telluric is None:
            self.flash_status_message("Bad value.")
            return

        self.active_spectrum.velocity_telluric = vel_telluric
        if resolution <= 0:
            resolution = None
        else:
            self.active_spectrum.resolution_telluric = resolution
            self.active_spectrum.resolution_atomic = resolution

        # Remove drawd lines if they exist
        self.remove_drawn_fitted_lines()


        self.operation_in_progress = True
        self.status_message("Fitting lines...")
        thread = threading.Thread(target=self.on_fit_lines_thread, args=(resolution, vel_telluric, selected_linelist, max_atomic_wave_diff, free_mu, check_derivatives))
        thread.setDaemon(True)
        thread.start()

    def on_fit_lines_thread(self, resolution, vel_telluric, selected_linelist, max_atomic_wave_diff, free_mu, check_derivatives):
        self.__update_numpy_arrays_from_widgets("lines")

        if resolution is not None and resolution > 0:
            logging.info("Smoothing spectrum...")
            self.queue.put((self.status_message, ["Smoothing spectrum..."], {}))
            smoothed_spectrum = ispec.convolve_spectrum(self.active_spectrum.data, 2*resolution, frame=self)
        else:
            smoothed_spectrum = None

        logging.info("Fitting lines...")
        self.queue.put((self.status_message, ["Fitting lines..."], {}))

        linemasks = ispec.fit_lines(self.regions["lines"], self.active_spectrum.data, self.active_spectrum.continuum_model, \
                            atomic_linelist=self.atomic_linelist[selected_linelist], \
                            max_atomic_wave_diff = max_atomic_wave_diff, \
                            telluric_linelist=self.telluric_linelist, vel_telluric=vel_telluric, \
                            discard_gaussian=False, discard_voigt=True, \
                            check_derivatives=check_derivatives, smoothed_spectrum=smoothed_spectrum, \
                            free_mu=free_mu, crossmatch_with_mu=free_mu, closest_match=False)
        # Exclude lines that have not been successfully cross matched with the atomic data
        # because we cannot calculate the chemical abundance (it will crash the corresponding routines)
        rejected_by_atomic_line_not_found = (linemasks['wave_nm'] == 0)
        linemasks = linemasks[~rejected_by_atomic_line_not_found]

        self.queue.put((self.on_fit_lines_finnish, [linemasks], {}))

    def on_fit_lines_finnish(self, linemasks, conserve_previous_regions=True):
        elements = "lines"

        if linemasks is not None and conserve_previous_regions and len(self.region_widgets["lines"]) - len(linemasks) > 0:
            diff_num_regions = len(self.region_widgets["lines"]) - len(linemasks)
            # Find regions that has been discarded due to a bad fit or other reason
            # in order to recover them
            recovered_regions = np.recarray((diff_num_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|U100')])
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
                    #recovered_regions[i]['note'] = region.get_note_text()
                    recovered_regions[i]['note'] = ""
                    i += 1
        else:
            diff_num_regions = 0

        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements, check_not_saved=False)

        self.active_spectrum.linemasks = linemasks

        if linemasks is not None and len(linemasks) > 0:
            total_regions = len(linemasks)
            line_regions = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|U100')])
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
                line_extra = str(line['wave_nm']) + ";" + str(line['spectrum_moog_species']) + ";"
                line_extra = line_extra + str(line['lower_state_cm1']) + ";" + str(line['upper_state_cm1']) + ";"
                line_extra = line_extra + str(line['loggf']) + ";" + str(line['spectrum_fudge_factor']) + ";"
                line_extra = line_extra + str(line['spectrum_transition_type']) + ";" + str(line['rad']) + ";"
                line_extra = line_extra + str(line['stark']) + ";" + str(line['waals']) + ";"
                line_extra = line_extra + str(line['ew']) + ";" + line['element'] + ";"
                line_extra = line_extra + str(line['telluric_wave_peak']) + ";" + str(line['telluric_depth'])
                line_extras.append(line_extra)
                line_model = ispec.GaussianModel(baseline=line['baseline'], A=line['A'], sig=line['sig'], mu=line['mu'])
                line_model.set_emu(line['mu_err'])
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
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = FindContinuumDialog(self, "Properties for finding continuum regions", self.find_continuum_regions_wave_step, self.find_continuum_regions_sigma, self.find_continuum_regions_max_continuum_diff)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        resolution = None
        fixed_wave_step = self.active_spectrum.dialog[key].results["Check for regions of minimum size"]
        sigma = self.active_spectrum.dialog[key].results["Maximum standard deviation"]
        max_continuum_diff = self.active_spectrum.dialog[key].results["Maximum fitted continuum difference (%)"]
        in_segments = self.active_spectrum.dialog[key].results["Look for continuum regions in"] == "Only inside segments"
        self.active_spectrum.dialog[key].destroy()

        if fixed_wave_step is None or sigma is None or max_continuum_diff is None:
            self.queue.put((self.flash_status_message, ["Bad value."], {}))
            return
        # Save values
        self.find_continuum_regions_wave_step = fixed_wave_step
        self.find_continuum_regions_sigma = sigma
        self.find_continuum_regions_max_continuum_diff = max_continuum_diff
        # Convert from % to over 1
        max_continuum_diff = max_continuum_diff/ 100

        if in_segments and (self.region_widgets["segments"] is None or len(self.region_widgets["segments"]) == 0):
            self.queue.put((self.flash_status_message, ["No segments found."], {}))
            return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_find_continuum_thread, args=(resolution, sigma, max_continuum_diff), kwargs={'fixed_wave_step':fixed_wave_step, 'in_segments':in_segments})
        thread.setDaemon(True)
        thread.start()

    def update_numpy_arrays_from_widgets(self, elements):
        total_regions = len(self.region_widgets[elements])
        if elements == "lines":
            self.regions[elements] = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|U100')])
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
        self.queue.put((self.status_message, ["Finding continuum regions..."], {}))
        if in_segments:
            self.update_numpy_arrays_from_widgets("segments")
            continuum_regions = ispec.find_continuum(self.active_spectrum.data, resolution, segments=self.regions["segments"], max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)
        else:
            continuum_regions = ispec.find_continuum(self.active_spectrum.data, resolution, max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)

        self.queue.put((self.on_find_continuum_finish, [continuum_regions], {}))

    def on_find_continuum_finish(self, continuum_regions):
        elements = "continuum"
        self.regions[elements] = continuum_regions
        self.draw_regions(elements)
        self.not_saved[elements] = True
        self.update_title()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Automatic finding of continuum regions ended.")

    def on_adjust_lines(self):
        if not self.check_active_spectrum_exists():
            return
        #if not self.check_continuum_model_exists():
            #return
        if self.check_operation_in_progress():
            return

        if self.not_saved["lines"]:
            msg = "Are you sure you want to modify the line masks without saving the current ones?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.update_numpy_arrays_from_widgets("lines")
        if self.regions["lines"] is None or len(self.regions["lines"]) <= 0:
            msg = 'There are no line masks'
            title = "Missing line masks"
            self.error(title, msg)
            self.flash_status_message("No line masks.")

        if self.active_spectrum.resolution_telluric == 0.0:
            R = self.active_spectrum.resolution_atomic
        else:
            R = self.active_spectrum.resolution_telluric

        key = "AdjustLinesDialog"
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = AdjustLinesDialog(self, "Adjusting line masks", margin=0.5, resolution=R)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        margin = self.active_spectrum.dialog[key].results["Margin around lines"]
        resolution = self.active_spectrum.dialog[key].results["Resolution"]
        check_derivatives = self.active_spectrum.dialog[key].results["Check derivatives before fitting"] == 1
        self.active_spectrum.dialog[key].destroy()

        if resolution <= 0:
            resolution = None
        else:
            self.active_spectrum.resolution_telluric = resolution
            self.active_spectrum.resolution_atomic = resolution

        elements = "lines"
        linemasks = self.regions["lines"]

        if resolution is not None:
            logging.info("Smoothing spectrum...")
            #self.queue.put((self.status_message, ["Smoothing spectrum..."], {}))
            smoothed_spectrum = ispec.convolve_spectrum(self.active_spectrum.data, 2*resolution, frame=self)
        else:
            smoothed_spectrum = self.active_spectrum.data

        logging.info("Adjusting line masks...")
        linemasks = ispec.adjust_linemasks(smoothed_spectrum, linemasks, max_margin=margin, check_derivatives=check_derivatives)

        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements, check_not_saved=False)
        self.regions[elements] = linemasks
        self.draw_regions(elements)
        self.not_saved[elements] = False
        self.update_title()
        self.flash_status_message("Line mask adjusted")
        self.canvas.draw()

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

        vel_telluric = self.active_spectrum.velocity_telluric


        key = "FindLinesDialog"
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = FindLinesDialog(self, "Properties for finding line masks", self.find_lines_min_depth, self.find_lines_max_depth, vel_telluric, R, "Fe 1, Fe 2", self.lists, self.default_lists)
        self.active_spectrum.dialog[key].show(updated_vel_telluric=vel_telluric)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        min_depth = self.active_spectrum.dialog[key].results["Minimum depth (% of the continuum)"]
        max_depth = self.active_spectrum.dialog[key].results["Maximum depth (% of the continuum)"]
        elements = self.active_spectrum.dialog[key].results["Select elements (comma separated)"]
        resolution = self.active_spectrum.dialog[key].results["Resolution"]
        vel_telluric = self.active_spectrum.dialog[key].results["Velocity respect to telluric lines (km/s)"]
        discard_tellurics = self.active_spectrum.dialog[key].results["Discard affected by tellurics"] == 1
        check_derivatives = self.active_spectrum.dialog[key].results["Check derivatives before fitting"] == 1
        in_segments = self.active_spectrum.dialog[key].results["Look for line masks in"] == "Only inside segments"
        selected_linelist = self.active_spectrum.dialog[key].results["Line list"]
        ### Find filename
        i = np.where(self.lists['atomic_lines']['name'] == selected_linelist)
        atomic_linelist_file = self.lists['atomic_lines']['path'][i][0]
        ####
        max_atomic_wave_diff = self.active_spectrum.dialog[key].results["Maximum atomic wavelength difference"]
        self.active_spectrum.dialog[key].destroy()

        telluric_linelist_file = resource_path("input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")

        # Read
        if self.molecules is None:
            self.molecules = ispec.read_molecular_symbols(molecules_file)
        if self.chemical_elements is None:
            self.chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
        if not selected_linelist in list(self.atomic_linelist.keys()):
            self.atomic_linelist[selected_linelist] = ispec.read_atomic_linelist(atomic_linelist_file)
            logging.warning("Limiting linelist to lines that have at least 0.01 depth in the Sun")
            solar = self.atomic_linelist[selected_linelist]['theoretical_depth'] >= 0.01
            self.atomic_linelist[selected_linelist] = self.atomic_linelist[selected_linelist][solar]
        if self.telluric_linelist is None:
            self.telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

        if max_depth is None or min_depth is None or resolution is None or vel_telluric is None or max_depth <= min_depth or max_depth <= 0 or min_depth < 0 or resolution < 0:
            self.queue.put((self.flash_status_message, ["Bad value."], {}))
            return
        ## Save values
        self.find_lines_max_depth = max_depth
        self.find_lines_min_depth = min_depth
        self.active_spectrum.velocity_telluric = vel_telluric
        if resolution <= 0:
            resolution = None
        else:
            self.active_spectrum.resolution_telluric = resolution
            self.active_spectrum.resolution_atomic = resolution

        if in_segments and (self.region_widgets["segments"] is None or len(self.region_widgets["segments"]) == 0):
            self.queue.put((self.flash_status_message, ["No segments found."], {}))
            return

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_find_lines_thread, args=(max_depth, min_depth, elements, resolution, vel_telluric, discard_tellurics, selected_linelist, max_atomic_wave_diff, check_derivatives), kwargs={'in_segments':in_segments})
        thread.setDaemon(True)
        thread.start()

    def on_find_lines_thread(self, max_depth, min_depth, elements, resolution, vel_telluric, discard_tellurics, selected_linelist, max_atomic_wave_diff, check_derivatives, in_segments=False):
        if in_segments:
            # Select spectrum from regions
            self.update_numpy_arrays_from_widgets("segments")
            spectrum = None
            for region in self.regions["segments"]:
                wave_base = region["wave_base"]
                wave_top = region["wave_top"]

                wfilter = np.logical_and(self.active_spectrum.data['waveobs'] >= wave_base, self.active_spectrum.data['waveobs'] <= wave_top)
                if spectrum is None:
                    spectrum = self.active_spectrum.data[wfilter]
                else:
                    spectrum = np.hstack((spectrum, self.active_spectrum.data[wfilter]))
        else:
            spectrum = self.active_spectrum.data

        if resolution is not None and resolution > 0:
            logging.info("Smoothing spectrum...")
            self.queue.put((self.status_message, ["Smoothing spectrum..."], {}))
            smoothed_spectrum = ispec.convolve_spectrum(spectrum, 2*resolution, frame=self)
        else:
            smoothed_spectrum = None

        self.queue.put((self.status_message, ["Generating line masks, fitting gaussians and matching atomic lines..."], {}))
        logging.info("Generating line masks, fitting gaussians and matching atomic lines...")

        linemasks = ispec.find_linemasks(spectrum, self.active_spectrum.continuum_model, \
                                atomic_linelist=self.atomic_linelist[selected_linelist], \
                                max_atomic_wave_diff = max_atomic_wave_diff, \
                                telluric_linelist=self.telluric_linelist, \
                                vel_telluric=vel_telluric, \
                                minimum_depth=min_depth, maximum_depth=max_depth, \
                                smoothed_spectrum=smoothed_spectrum, \
                                check_derivatives=check_derivatives, \
                                discard_gaussian=False, discard_voigt=True, \
                                closest_match=False, \
                                frame=self)

        # If no peaks found, just finnish
        if linemasks is None or len(linemasks) == 0:
            self.queue.put((self.on_find_lines_finish, [None], {}))
            return

        logging.info("Applying filters to discard bad line masks...")
        self.queue.put((self.status_message, ["Applying filters to discard bad line masks..."], {}))

        rejected_by_atomic_line_not_found = (linemasks['wave_nm'] == 0)
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
            # Check that there are no linemasks that exceed segment limits (because
            # the spectrum has been cutted, the algorithm can create this overdimensioned linemasks)
            # If so, change the linemask limit to the segment limit.
            for region in self.regions["segments"]:
                wave_base = region["wave_base"]
                wave_top = region["wave_top"]
                wfilter = np.logical_and(linemasks['wave_peak'] >= wave_base, linemasks['wave_peak'] <= wave_top)
                filter_base_exceeded = np.where(np.logical_and(wfilter, linemasks['wave_base'] < wave_base))[0]
                filter_top_exceeded = np.where(np.logical_and(wfilter, linemasks['wave_top'] > wave_top))[0]
                if len(filter_base_exceeded) > 0:
                    linemasks['wave_base'][filter_base_exceeded] = wave_base
                if len(filter_top_exceeded) > 0:
                    linemasks['wave_top'][filter_top_exceeded] = wave_top

        linemasks = linemasks[~discarded]
        total_regions = len(linemasks)

        # If no regions found, just finnish
        if total_regions == 0:
            self.queue.put((self.on_find_lines_finish, [None], {}))
            return


        self.queue.put((self.on_find_lines_finish, [linemasks], {}))

    def on_find_lines_finish(self, linemasks):
        self.on_fit_lines_finnish(linemasks, conserve_previous_regions=False)


    def on_create_segments_around_lines(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        self.update_numpy_arrays_from_widgets("lines")
        if self.regions["lines"] is None or len(self.regions["lines"]) <= 0:
            msg = 'There are no line masks'
            title = "Missing line masks"
            self.error(title, msg)
            self.flash_status_message("No line masks.")
            return

        self.update_numpy_arrays_from_widgets("segments")
        if self.regions["segments"] is not None and \
                len(self.regions["segments"]) > 0 and \
                self.not_saved["segments"] == True:
            msg = "Are you sure you want to continue without saving the current segments?"
            title = "Segments not saved"
            if not self.question(title, msg):
                return

        key = "FindSegmentsDialog"
        if key not in self.dialog:
            self.dialog[key] = FindSegmentsDialog(self, "Properties for finding segments", margin=0.25)
        self.dialog[key].show()

        if self.dialog[key].results is None:
            self.dialog[key].destroy()
            return

        margin = self.dialog[key].results["Margin around lines"]
        self.dialog[key].destroy()


        elements = "segments"
        self.remove_regions(elements)

        linemasks = self.regions["lines"]
        self.regions[elements] = ispec.create_segments_around_lines(linemasks, margin=margin)

        self.draw_regions(elements)
        self.not_saved[elements] = False
        self.update_title()
        self.flash_status_message("Segments around line mask created")
        self.canvas.draw()


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
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = CorrectVelocityDialog(self, "Velocity correction", vel_type, default_vel)
        self.active_spectrum.dialog[key].show(updated_vel=default_vel)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return


        velocity = self.active_spectrum.dialog[key].results["Velocity relative to %s (km/s)" % vel_type]
        in_regions = self.active_spectrum.dialog[key].results["Apply correction on"] == "Regions"
        self.active_spectrum.dialog[key].destroy()

        if velocity is None:
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
                        self.regions[elements] = ispec.correct_velocity_regions(self.regions[elements], velocity, with_peak=True)
                    else:
                        self.regions[elements] = ispec.correct_velocity_regions(self.regions[elements], velocity)
                    self.draw_regions(elements)
                    self.not_saved[elements] = False
        else:
            # Modify velocities according to the correction
            self.active_spectrum.velocity_atomic -= velocity
            self.active_spectrum.velocity_telluric -= velocity
            self.active_spectrum.velocity_template -= velocity
            # Correct
            self.active_spectrum.data = ispec.correct_velocity(self.active_spectrum.data, velocity)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a " + vel_type + " velocity correction of %s." % velocity)


    def on_determine_barycentric_vel(self):
        key = "DetermineBarycentricCorrectionDialog"
        if key not in self.dialog:
            self.dialog[key] = DetermineBarycentricCorrectionDialog(self,  "Barycentric velocity determination", "15/02/2012", "00:00:00", "19:50:46.99", "08:52:5.96")
        self.dialog[key].show()

        if self.dialog[key].results is None:
            # Cancel
            self.dialog[key].destroy()
            return

        date_string = self.dialog[key].results["Date (DD/MM/YYY)"]
        time_string = self.dialog[key].results["Time (HH:MM:SS)"]
        ra = self.dialog[key].results["Right ascension (HH:MM:SS)"]
        dec = self.dialog[key].results["Declination (DD:MM:SS)"]
        self.dialog[key].destroy()

        try:
            day, month, year = list(map(float, date_string.split("/")))
            hours, minutes, seconds = list(map(float, time_string.split(":")))
            ra_hours, ra_minutes, ra_seconds = list(map(float, ra.split(":")))
            dec_degree, dec_minutes, dec_seconds = list(map(float, dec.split(":")))
        except:
            msg = 'Some input values are not in the expected data format.'
            title = "Bad values"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        if None in [day, month, year, hours, minutes, seconds, ra_hours, ra_minutes, ra_seconds, dec_degree, dec_minutes, dec_seconds]:
            self.flash_status_message("Bad value.")
            return

        # Project velocity toward star
        self.barycentric_vel = ispec.calculate_barycentric_velocity_correction((year, month, day, hours, minutes, seconds), (ra_hours, ra_minutes, ra_seconds, dec_degree, dec_minutes, dec_seconds))

        msg = "Barycentric velocity determined: " + str(self.barycentric_vel) + " km/s"
        title = "Barycentric velocity"
        self.info(title, msg)
        self.flash_status_message(msg)

    def on_estimate_snr(self):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.snr is not None:
            msg = "Previous estimated SNR: %.2f. Re-estimate again? " % self.active_spectrum.snr
            title = "Signal-to-Noise Ratio"
            if not self.question(title, msg):
                return

        key = "EstimateSNRDialog"
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = EstimateSNRDialog(self, "Properties for estimating SNR", num_points=10)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        num_points = self.active_spectrum.dialog[key].results["Number of points"]
        wave_step = self.active_spectrum.dialog[key].results["Wavelength step (resampling)"]
        estimate_from_flux = self.active_spectrum.dialog[key].results["Estimate SNR"] == "From fluxes in blocks of N points"
        self.active_spectrum.dialog[key].destroy()

        if num_points is None or wave_step is None:
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
                estimated_snr = np.median(spec['flux']/ spec['err'])
                self.on_estimate_snr_finnish(estimated_snr)
            else:
                msg = 'All value errors are set to zero or negative numbers'
                title = "SNR estimation error"
                self.error(title, msg)


    def on_estimate_snr_thread(self, num_points, wave_step):
        self.queue.put((self.status_message, ["Resampling spectrum..."], {}))
        self.queue.put((self.update_progress, [10], {}))
        xaxis = np.arange(np.min(self.active_spectrum.data["waveobs"]), np.max(self.active_spectrum.data["waveobs"]), wave_step)
        resampled_spectrum_data = ispec.resample_spectrum(self.active_spectrum.data, xaxis, frame=self)

        self.queue.put((self.status_message, ["Estimating SNR for the whole spectrum..."], {}))
        estimated_snr = ispec.estimate_snr(resampled_spectrum_data['flux'], num_points=num_points, frame=self)
        self.queue.put((self.on_estimate_snr_finnish, [estimated_snr], {}))

    def on_estimate_snr_finnish(self, estimated_snr):
        self.active_spectrum.snr = estimated_snr
        msg = "Estimated SNR: %.2f" % self.active_spectrum.snr
        title = "Signal-to-Noise Ratio"
        self.info(title, msg)
        self.flash_status_message(msg)
        self.operation_in_progress = False

    def on_calculate_errors(self):
        if self.check_operation_in_progress():
            return

        if np.all(self.active_spectrum.data['err'] > 0.0) and np.all(self.active_spectrum.data['err'] < np.max(self.active_spectrum.data['flux'])):
            mean_err = np.mean(self.active_spectrum.data['err'])
            max_err = np.max(self.active_spectrum.data['err'])
            min_err = np.min(self.active_spectrum.data['err'])
            msg = "Spectrum seems to already have errors (mean: %.4f, max: %.4f, min: %.4f). Are you sure to overwrite them with new estimates? " % (mean_err, max_err, min_err)
            title = "Spectrum's errors"
            if not self.question(title, msg):
                return

        key = "EstimateErrorsDialog"
        if key not in self.active_spectrum.dialog:
            snr = 10.0
            self.active_spectrum.dialog[key] = EstimateErrorsDialog(self, "Calculate spectrum errors", snr)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        snr = self.active_spectrum.dialog[key].results["SNR (Signal-to-Noise Ratio)"]
        self.active_spectrum.dialog[key].destroy()

        if snr is None or snr <= 0:
            msg = "SNR should be greater than zero"
            title = "SNR error"
            self.error(title, msg)
            return

        self.active_spectrum.data['err'] = self.active_spectrum.data['flux']/ snr
        self.active_spectrum.not_saved = True

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        msg = "Spectrum's errors estimated!"
        self.flash_status_message(msg)

    def on_add_noise(self):
        if self.check_operation_in_progress():
            return

        key = "AddNoiseDialog"
        if key not in self.active_spectrum.dialog:
            snr = 10.0
            self.active_spectrum.dialog[key] = AddNoiseDialog(self, "Add noise to spectrum", snr)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        snr = self.active_spectrum.dialog[key].results["SNR (Signal-to-Noise Ratio)"]
        distribution = self.active_spectrum.dialog[key].results["Distribution"]
        self.active_spectrum.dialog[key].destroy()

        if snr is None or snr <= 0:
            msg = "SNR should be greater than zero"
            title = "SNR error"
            self.error(title, msg)
            return

        self.active_spectrum.data = ispec.add_noise(self.active_spectrum.data, snr, distribution)
        self.active_spectrum.not_saved = True

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        msg = "Noise added to the spectrum!"
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
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = DegradeResolutionDialog(self, "Degrade spectrum resolution", R, R/2)
        self.active_spectrum.dialog[key].show(updated_from_resolution=R, updated_to_resolution=R/2)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        from_resolution = self.active_spectrum.dialog[key].results["Initial resolution"]
        to_resolution = self.active_spectrum.dialog[key].results["Final resolution"]
        self.active_spectrum.dialog[key].destroy()

        if from_resolution is None or to_resolution is None or (from_resolution != 0 and from_resolution <= to_resolution) or from_resolution < 0 or to_resolution <= 0:
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
        self.queue.put((self.status_message, ["Degrading spectrum resolution..."], {}))
        if from_resolution == 0:
            # Smooth
            convolved_spectrum = ispec.convolve_spectrum(self.active_spectrum.data, to_resolution, frame=self)
        else:
            convolved_spectrum = ispec.convolve_spectrum(self.active_spectrum.data, to_resolution, from_resolution=from_resolution, frame=self)
        self.queue.put((self.on_degrade_resolution_finnish, [convolved_spectrum, from_resolution, to_resolution], {}))

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
        if key not in self.dialog:
            self.dialog[key] = CleanSpectrumDialog(self, "Clean fluxes and errors", 0.0, 1.2, 0.0, 1.2)
        self.dialog[key].show()

        if self.dialog[key].results is None:
            # Cancel
            self.dialog[key].destroy()
            return

        filter_by_flux = self.dialog[key].results["Filter by flux"] == 1
        flux_base = self.dialog[key].results["Base flux"]
        flux_top = self.dialog[key].results["Top flux"]
        filter_by_error = self.dialog[key].results["Filter by error"] == 1
        err_base = self.dialog[key].results["Base error"]
        err_top = self.dialog[key].results["Top error"]
        filter_cosmics = self.dialog[key].results["Filter cosmics"] == 1
        resampling_step = self.dialog[key].results["Resampling step"]
        window_size = self.dialog[key].results["Window size"]
        variation_limit = self.dialog[key].results["Variation limit"]

        replace_by = self.dialog[key].results["Replace by"]
        self.dialog[key].destroy()

        if (filter_cosmics or replace_by == "Continuum") and self.active_spectrum.continuum_model is None:
            msg = "It is necessary to previously fit continuum in order to clean values with the selected options"
            title = "Continuum not fitted"
            self.error(title, msg)
            self.flash_status_message("Continuum not fitted.")
            return

        if flux_base is None or flux_top is None or flux_top <= flux_base \
                or err_base is None or err_top is None or err_top <= err_base:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to clean it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Cleaning spectrum...")

        if filter_by_flux and filter_by_error and filter_cosmics:
            cosmics = ispec.create_filter_cosmic_rays(self.active_spectrum.data, self.active_spectrum.continuum_model, resampling_wave_step=resampling_step, window_size=window_size, variation_limit=variation_limit)
            cfilter = np.logical_not(cosmics)
            ffilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
            efilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
            wfilter = np.logical_and(ffilter, efilter)
            wfilter = np.logical_and(wfilter, cfilter)
        if filter_by_flux and filter_by_error:
            ffilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
            efilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
            wfilter = np.logical_and(ffilter, efilter)
        elif filter_by_flux:
            wfilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
        elif filter_by_error:
            wfilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
        elif filter_cosmics:
            cosmics = ispec.create_filter_cosmic_rays(self.active_spectrum.data, self.active_spectrum.continuum_model, resampling_wave_step=resampling_step, window_size=window_size, variation_limit=variation_limit)
            wfilter = np.logical_not(cosmics)
        else:
            wfilter = np.logical_not(np.isnan(self.active_spectrum.data['err']))

        if len(self.active_spectrum.data[wfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong flux/error ranges"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        if len(self.active_spectrum.data[wfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong ranges"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return

        if replace_by == "Zeros":
            self.active_spectrum.data['flux'][~wfilter] = 0.0
            self.active_spectrum.data['err'][~wfilter] = 0.0
        elif replace_by == "NaN":
            self.active_spectrum.data['flux'][~wfilter] = np.nan
            self.active_spectrum.data['err'][~wfilter] = np.nan
        elif replace_by == "Continuum":
            self.active_spectrum.data['flux'][~wfilter] = self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'][~wfilter])
            self.active_spectrum.data['err'][~wfilter] = 0.0
        else:
            self.active_spectrum.data = self.active_spectrum.data[wfilter]

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

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
        vel_telluric = self.active_spectrum.velocity_telluric
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = CleanTelluricsDialog(self, "Clean telluric regions", vel_telluric, -30.0, 30.0, 0.02)
        self.active_spectrum.dialog[key].show(updated_vel=vel_telluric)

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return

        vel_telluric = self.active_spectrum.dialog[key].results["Velocity relative to tellurics"]
        min_vel = self.active_spectrum.dialog[key].results["Minimum velocity"]
        max_vel = self.active_spectrum.dialog[key].results["Maximum velocity"]
        min_depth = self.active_spectrum.dialog[key].results["Minimum tellurics depth"]
        replace_by = self.active_spectrum.dialog[key].results["Replace by"]
        self.active_spectrum.dialog[key].destroy()

        if replace_by == "Continuum" and self.active_spectrum.continuum_model is None:
            msg = "It is necessary to previously fit continuum in order to clean values and replace them by continuum"
            title = "Continuum not fitted"
            self.error(title, msg)
            self.flash_status_message("Continuum not fitted.")
            return

        if vel_telluric is None or min_depth is None or min_vel is None or max_vel is None or min_vel >= max_vel:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to clean it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        self.status_message("Cleaning spectrum...")

        self.active_spectrum.velocity_telluric = vel_telluric

        if self.telluric_linelist is None:
            telluric_linelist_file = resource_path("input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst")
            self.telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.01)

        # - Filter regions that may be affected by telluric lines
        #dfilter = self.telluric_linelist['depth'] > np.percentile(self.telluric_linelist['depth'], 75) # (only the 25% of the deepest ones)
        dfilter = self.telluric_linelist['depth'] > min_depth
        tfilter = ispec.create_filter_for_regions_affected_by_tellurics(self.active_spectrum.data['waveobs'], \
                                    self.telluric_linelist[dfilter], min_velocity=-vel_telluric+min_vel, max_velocity=-vel_telluric+max_vel)

        if len(self.active_spectrum.data[tfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong ranges"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
            return
        if replace_by == "Zeros":
            self.active_spectrum.data['flux'][tfilter] = 0.0
            self.active_spectrum.data['err'][tfilter] = 0.0
        elif replace_by == "NaN":
            self.active_spectrum.data['flux'][tfilter] = np.nan
            self.active_spectrum.data['err'][tfilter] = np.nan
        elif replace_by == "Continuum":
            self.active_spectrum.data['flux'][tfilter] = self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'][tfilter])
            self.active_spectrum.data['err'][tfilter] = 0.0
        else:
            self.active_spectrum.data = self.active_spectrum.data[~tfilter]

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

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
        if key not in self.dialog:
            self.dialog[key] = CutSpectrumDialog(self, "Wavelength range reduction", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2))
        self.dialog[key].show()

        if self.dialog[key].results is None:
            # Cancel
            self.dialog[key].destroy()
            return

        wave_base = self.dialog[key].results["Base wavelength"]
        wave_top = self.dialog[key].results["Top wavelength"]
        in_segments = self.dialog[key].results["Consider"] == "Segments"
        in_linemasks = self.dialog[key].results["Consider"] == "Line masks"
        replace_by = self.dialog[key].results["Replace by"]
        self.dialog[key].destroy()

        if (not in_segments and not in_linemasks) and (wave_base is None or wave_top is None or wave_top <= wave_base):
            self.flash_status_message("Bad value.")
            return

        if in_segments and (self.region_widgets["segments"] is None or len(self.region_widgets["segments"]) == 0):
            self.queue.put((self.flash_status_message, ["No segments found."], {}))
            return

        if in_linemasks and (self.region_widgets["lines"] is None or len(self.region_widgets["lines"]) == 0):
            self.queue.put((self.flash_status_message, ["No line masks found."], {}))
            return

        if not in_segments and not in_linemasks:
            wfilter = ispec.create_wavelength_filter(self.active_spectrum.data, wave_base=wave_base, wave_top=wave_top)
        elif in_segments:
            self.__update_numpy_arrays_from_widgets("segments")
            wfilter = ispec.create_wavelength_filter(self.active_spectrum.data, regions=self.regions["segments"])
        else:
            self.__update_numpy_arrays_from_widgets("lines")
            wfilter = ispec.create_wavelength_filter(self.active_spectrum.data, regions=self.regions["lines"])

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

        self.active_spectrum.data['flux'][~wfilter] = 0.0
        self.active_spectrum.data['err'][~wfilter] = 0.0
        if replace_by == "Zeros":
            self.active_spectrum.data['flux'][~wfilter] = 0.0
            self.active_spectrum.data['err'][~wfilter] = 0.0
        elif replace_by == "NaN":
            self.active_spectrum.data['flux'][~wfilter] = np.nan
            self.active_spectrum.data['err'][~wfilter] = np.nan
        elif replace_by == "Continuum":
            self.active_spectrum.data['flux'][~wfilter] = self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'][~wfilter])
            self.active_spectrum.data['err'][~wfilter] = 0.0
        else:
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
        if key not in self.dialog:
            self.dialog[key] = ResampleSpectrumDialog(self, "Resample spectrum", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2), 0.001, median_step, mean_step, min_step, max_step)
        self.dialog[key].show(updated_median_step=median_step, updated_mean_step=mean_step, updated_min_step=min_step, updated_max_step=max_step)

        if self.dialog[key].results is None:
            # Cancel
            self.dialog[key].destroy()
            return

        wave_base = self.dialog[key].results["Base wavelength"]
        wave_top = self.dialog[key].results["Top wavelength"]
        wave_step = self.dialog[key].results["Wavelength step"]
        method = self.dialog[key].results["Method"]
        self.dialog[key].destroy()

        if wave_base is None or wave_top is None or wave_top <= wave_base or wave_step <= 0:
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
        self.queue.put((self.status_message, ["Resampling spectrum..."], {}))
        self.queue.put((self.update_progress, [10], {}))
        xaxis = np.arange(wave_base, wave_top, wave_step)
        resampled_spectrum_data = ispec.resample_spectrum(self.active_spectrum.data, xaxis, method=method, frame=self)
        self.active_spectrum.data = resampled_spectrum_data
        self.queue.put((self.on_resample_spectrum_finnish, [], {}))

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
            if wave_base is None:
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

        if dialog.results is None:
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


        if wave_base is None or wave_top is None or wave_top <= wave_base or wave_step <= 0:
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
            self.queue.put((self.status_message, ["Resampling spectrum %i of %i (%s)..." % (i+1, total, spec.name)], {}))
            if spec == self.active_spectrum:
                active = i
            resampled_spectrum_data = ispec.resample_spectrum(spec.data, xaxis, frame=self)
            resampled_spectra.append(resampled_spectrum_data)
            i += 1

        # Combine
        self.queue.put((self.status_message, ["Combining spectra..."], {}))
        self.queue.put((self.update_progress, [10], {}))

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
                combined_spectrum = ispec.create_spectrum_structure(xaxis, err=std)
                combined_spectrum['flux'] = np.mean(matrix, axis=0)
                combined_spectrum_name = "Cumulative_mean_spectrum"
            else:
                # Median fluxes
                combined_spectrum = ispec.create_spectrum_structure(xaxis, err=std)
                combined_spectrum['flux'] = np.median(matrix, axis=0)
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
            combined_spectrum = ispec.create_spectrum_structure(xaxis, flux, err)
            combined_spectrum_name = "Subtracted_spectrum"
        elif operation_add:
            flux = np.zeros(total_wavelengths)
            err = np.zeros(total_wavelengths)
            for spec in resampled_spectra:
                flux = flux + spec['flux']
                # Error propagation assuming that they are independent
                err = np.sqrt(np.power(err,2) + np.power(spec['err'],2))
            combined_spectrum = ispec.create_spectrum_structure(xaxis, flux, err)
            combined_spectrum_name = "Added_spectrum"
        elif operation_divide:
            flux = resampled_spectra[active]['flux'].copy()
            err = resampled_spectra[active]['err'].copy()
            i = 0
            for spec in resampled_spectra:
                #logging.warning("Division by zero may occur")
                if i != active:
                    # Error propagation assuming that they are independent
                    err = np.sqrt(np.power(flux/ err, 2) + np.power(spec['flux']/ spec['err'], 2)) * 0.7
                    flux = flux * (1. / spec['flux'])
                i += 1
            combined_spectrum = ispec.create_spectrum_structure(xaxis, flux, err)
            combined_spectrum_name = "Divided_spectrum"

        # Free memory
        for spec in resampled_spectra:
            del spec
        del resampled_spectra

        self.queue.put((self.on_combine_spectra_finnish, [combined_spectrum, combined_spectrum_name], {}))

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
        #self.active_spectrum.data['flux'] /= self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'])
        #self.active_spectrum.data['err'] /= self.active_spectrum.continuum_model(self.active_spectrum.data['waveobs'])
        self.active_spectrum.data = ispec.normalize_spectrum(self.active_spectrum.data, self.active_spectrum.continuum_model)
        self.active_spectrum.not_saved = True

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Establish the new continuum at 1.0
        self.active_spectrum.continuum_model = ispec.fit_continuum(self.active_spectrum.data, fixed_value=1.0, model="Fixed value")
        waveobs = self.active_spectrum.data['waveobs']
        self.active_spectrum.continuum_data = ispec.create_spectrum_structure(waveobs, self.active_spectrum.continuum_model(waveobs))
        self.draw_continuum_spectrum()

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum normalized!")


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
        if key not in self.dialog:
            operation_waveobs = self.operation_waveobs
            operation_flux = self.operation_flux
            operation_err = self.operation_err
            self.dialog[key] = OperateSpectrumDialog(self, "Apply mathematical expression", self.safe_operations_description, operation_waveobs, operation_flux, operation_err)
        self.dialog[key].show()

        if self.dialog[key].results is None:
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

    def on_interpolate(self):
        if self.check_operation_in_progress():
            return

        if len(self.lists['grid']) > 0:

            if self.active_spectrum is not None:
                wave_base = np.round(np.min(self.active_spectrum.data['waveobs']), 2)
                wave_top = np.round(np.max(self.active_spectrum.data['waveobs']), 2)
            else:
                wave_base = 515.0 # Magnesium triplet region
                wave_top = 525.0
                #wave_top = 517.0
            teff = 5771.0
            logg = 4.44
            MH = 0.00
            alpha = 0.00
            macroturbulence = 4.21
            vsini = 1.6
            limb_darkening_coeff = 0.6
            microturbulence_vel = 1.07
            #resolution = 47000
            #resolution = 0
            resolution = 300000
            wave_step = 0.001

            key = "InterpolateSpectrumDialog"
            if key not in self.dialog:
                self.dialog[key] = InterpolateSpectrumDialog(self, "Spectrum interpolator", wave_base, wave_top, wave_step, resolution, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, self.lists, self.default_lists)
            self.dialog[key].show()

            if self.dialog[key].results is None:
                # Cancel
                self.dialog[key].destroy()
                return

            selected_grid = self.dialog[key].results["Grid"]
            teff = self.dialog[key].results["Effective temperature (K)"]
            logg = self.dialog[key].results["Surface gravity (log g)"]
            MH = self.dialog[key].results["Metallicity [M/H]"]
            alpha = self.dialog[key].results["Alpha enhancement [alpha/Fe]"]
            microturbulence_vel = self.dialog[key].results["Microturbulence velocity (km/s)"]
            macroturbulence = self.dialog[key].results["Macroturbulence velocity (km/s)"]
            vsini = self.dialog[key].results["Rotation (v sin(i)) (km/s)"]
            limb_darkening_coeff = self.dialog[key].results["Limb darkening coefficient"]
            resolution = self.dialog[key].results["Resolution"]
            wave_base = self.dialog[key].results["Wavelength min (nm)"]
            wave_top = self.dialog[key].results["Wavelength max (nm)"]
            wave_step = self.dialog[key].results["Wavelength step (nm)"]
            in_segments = self.dialog[key].results["Generate spectrum for"] == "Segments"
            in_lines = self.dialog[key].results["Generate spectrum for"] == "Line masks"

            self.dialog[key].destroy()

            ### Find filenames
            i = np.where(self.lists['grid']['name'] == selected_grid)
            grid_dirname = os.path.dirname(self.lists['grid']['path'][i][0])

            if in_segments:
                elements_type = "segments"
            if in_lines:
                elements_type = "lines"

            if teff is None or logg is None or MH is None or alpha is None or microturbulence_vel is None or resolution is None or wave_base is None or wave_top is None:
                self.flash_status_message("Bad value.")
                return

            if selected_grid not in self.grid:
                logging.info("Loading grid %s..." % selected_grid)
                self.status_message("Loading grid %s..." % selected_grid)
                self.grid[selected_grid] = ispec.load_spectral_grid(grid_dirname)

            if not ispec.valid_interpolated_spectrum_target(self.grid[selected_grid], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic':microturbulence_vel}):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the grid limits."
                title = 'Out of the grid limits'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            if not in_segments and not in_lines:
                regions = None # Compute fluxes for all the wavelengths
            else:
                #in_segments or in_lines
                if len(self.region_widgets[elements_type]) == 0:
                    self.flash_status_message("No segments/line masks present for synthetic spectrum generation.")
                    return
                self.__update_numpy_arrays_from_widgets(elements_type)
                regions = self.regions[elements_type]
                wave_base = np.min(regions['wave_base'])
                wave_top = np.max(regions['wave_top'])

            if wave_base >= wave_top:
                msg = "Bad wavelength range definition, maximum value cannot be lower than minimum value."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            waveobs = np.arange(wave_base, wave_top, wave_step)
            if in_segments or in_lines:
                wfilter = ispec.create_wavelength_filter({'waveobs': waveobs}, wave_base=wave_base, wave_top=wave_top, regions=regions)
                waveobs = waveobs[wfilter]
            total_points = len(waveobs)

            if total_points < 2:
                msg = "Wavelength range too narrow."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            self.operation_in_progress = True
            self.status_message("Interpolating spectrum...")
            code = "grid"
            thread = threading.Thread(target=self.on_interpolate_thread, args=(code, selected_grid, waveobs, regions, teff, logg, MH, alpha, microturbulence_vel,  macroturbulence, vsini, limb_darkening_coeff, resolution, ))
            thread.setDaemon(True)
            thread.start()

    def on_interpolate_thread(self, code, selected_grid, waveobs, regions, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, resolution):

        synth_spectrum = ispec.create_spectrum_structure(waveobs)

        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])

        error_message = None
        try:
            # waveobs is multiplied by 10.0 in order to be converted from nm to armstrongs
            atmosphere_layers = None
            linelist = None
            isotopes = None
            abundances = None
            synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, R=resolution, regions=regions, verbose=1, gui_queue=self.queue, code=code, grid=self.grid[selected_grid])
        except Exception as e:
            error_message = str(e)


        synth_spectrum.sort(order='waveobs') # Make sure it is ordered by wavelength

        # Remove atmosphere model temporary file
        self.queue.put((self.on_synthesize_finnish, [code, synth_spectrum, teff, logg, MH, alpha, microturbulence_vel, error_message], {}))

    def on_synthesize_finnish(self, code, synth_spectrum, teff, logg, MH, alpha, microturbulence_vel, error_message):
        if error_message is not None:
            msg = error_message
            title = 'Problem synthesizing spectrum'
            self.error(title, msg)
            self.operation_in_progress = False
            return
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
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        # Name: If it already exists, add a suffix
        name = self.get_name(code + "_" + str(teff) + "_" + str(logg) + "_"  + str(MH) + "_"  + str(alpha) + "_" + str(microturbulence_vel))
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
        self.flash_status_message("Spectrum generated!")


    def on_synthesize(self):
        if self.check_operation_in_progress():
            return

        if (ispec.is_spectrum_support_enabled() \
                #or ispec.is_turbospectrum_support_enabled() \
                #or ispec.is_sme_support_enabled() \
                #or ispec.is_moog_support_enabled() \
                #or ispec.is_synthe_support_enabled() \
                ) and \
                len(self.lists['atmospheres']) > 0 and len(self.lists['abundances']) > 0 and len(self.lists['atomic_lines']) > 0:

            if self.active_spectrum is not None:
                wave_base = np.round(np.min(self.active_spectrum.data['waveobs']), 2)
                wave_top = np.round(np.max(self.active_spectrum.data['waveobs']), 2)
            else:
                wave_base = 515.0 # Magnesium triplet region
                wave_top = 525.0
                #wave_top = 517.0
            teff = 5771.0
            logg = 4.44
            MH = 0.00
            alpha = 0.00
            macroturbulence = 4.21
            vsini = 1.6
            limb_darkening_coeff = 0.6
            microturbulence_vel = 1.07
            #resolution = 47000
            #resolution = 0
            resolution = 300000
            wave_step = 0.001

            key = "SyntheticSpectrumDialog"
            if key not in self.dialog:
                self.dialog[key] = SyntheticSpectrumDialog(self, "Synthetic spectrum generator", wave_base, wave_top, wave_step, resolution, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, self.lists, self.default_lists)
            self.dialog[key].show()

            if self.dialog[key].results is None:
                # Cancel
                self.dialog[key].destroy()
                return

            code = self.dialog[key].results["Code"].lower()
            teff = self.dialog[key].results["Effective temperature (K)"]
            logg = self.dialog[key].results["Surface gravity (log g)"]
            MH = self.dialog[key].results["Metallicity [M/H]"]
            alpha = self.dialog[key].results["Alpha enhancement [alpha/Fe]"]
            microturbulence_vel = self.dialog[key].results["Microturbulence velocity (km/s)"]
            macroturbulence = self.dialog[key].results["Macroturbulence velocity (km/s)"]
            vsini = self.dialog[key].results["Rotation (v sin(i)) (km/s)"]
            limb_darkening_coeff = self.dialog[key].results["Limb darkening coefficient"]
            resolution = self.dialog[key].results["Resolution"]
            wave_base = self.dialog[key].results["Wavelength min (nm)"]
            wave_top = self.dialog[key].results["Wavelength max (nm)"]
            wave_step = self.dialog[key].results["Wavelength step (nm)"]
            selected_atmosphere_models = self.dialog[key].results["Model atmosphere"]
            selected_abundances = self.dialog[key].results["Solar abundances"]
            selected_linelist = self.dialog[key].results["Line list"]
            in_segments = self.dialog[key].results["Generate spectrum for"] == "Segments"
            in_lines = self.dialog[key].results["Generate spectrum for"] == "Line masks"

            self.dialog[key].destroy()

            ### Find filenames
            i = np.where(self.lists['atomic_lines']['name'] == selected_linelist)
            atomic_linelist_file = self.lists['atomic_lines']['path'][i][0]

            i = np.where(self.lists['abundances']['name'] == selected_abundances)
            abundances_file = self.lists['abundances']['path'][i][0]

            i = np.where(self.lists['atmospheres']['name'] == selected_atmosphere_models)
            atmospheres_file = os.path.dirname(self.lists['atmospheres']['path'][i][0])
            ####

            isotope_file = resource_path("input/isotopes/SPECTRUM.lst")

            if in_segments:
                elements_type = "segments"
            if in_lines:
                elements_type = "lines"

            if teff is None or logg is None or MH is None or alpha is None or microturbulence_vel is None or resolution is None or wave_base is None or wave_top is None:
                self.flash_status_message("Bad value.")
                return

            if selected_atmosphere_models not in self.modeled_layers_pack:
                logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
                self.modeled_layers_pack[selected_atmosphere_models] = ispec.load_modeled_layers_pack(atmospheres_file)

            if not ispec.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic':microturbulence_vel}):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models."
                title = 'Out of the atmospheric models'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            # Load SPECTRUM linelist
            chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
            molecules_file = resource_path("input/abundances/molecular_symbols.dat")
            if self.molecules is None:
                self.molecules = ispec.read_molecular_symbols(molecules_file)
            if self.chemical_elements is None:
                self.chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
            if not selected_linelist in list(self.atomic_linelist.keys()):
                self.atomic_linelist[selected_linelist] = ispec.read_atomic_linelist(atomic_linelist_file)
                logging.warning("Limiting linelist to lines that have at least 0.01 depth in the Sun")
                solar = self.atomic_linelist[selected_linelist]['theoretical_depth'] >= 0.01
                self.atomic_linelist[selected_linelist] = self.atomic_linelist[selected_linelist][solar]
            isotopes = ispec.read_isotope_data(isotope_file)
            linelist = self.atomic_linelist[selected_linelist]

            # Load SPECTRUM abundances
            if not abundances_file in list(self.solar_abundances.keys()):
                self.solar_abundances[abundances_file] = ispec.read_solar_abundances(abundances_file)
            abundances = self.solar_abundances[abundances_file]

            # Prepare atmosphere model
            self.status_message("Interpolating atmosphere model...")
            atmosphere_layers = ispec.interpolate_atmosphere_layers(self.modeled_layers_pack[selected_atmosphere_models], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha})


            if not in_segments and not in_lines:
                regions = None # Compute fluxes for all the wavelengths
            else:
                #in_segments or in_lines
                if len(self.region_widgets[elements_type]) == 0:
                    self.flash_status_message("No segments/line masks present for synthetic spectrum generation.")
                    return
                self.__update_numpy_arrays_from_widgets(elements_type)
                regions = self.regions[elements_type]
                wave_base = np.min(regions['wave_base'])
                wave_top = np.max(regions['wave_top'])

            if wave_base >= wave_top:
                msg = "Bad wavelength range definition, maximum value cannot be lower than minimum value."
                title = 'Wavelength range'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return
            waveobs = np.arange(wave_base, wave_top, wave_step)
            if in_segments or in_lines:
                wfilter = ispec.create_wavelength_filter({'waveobs': waveobs}, wave_base=wave_base, wave_top=wave_top, regions=regions)
                waveobs = waveobs[wfilter]

            # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectrum
            if np.min(waveobs) < 100.0 or np.max(waveobs) > 4000.0:
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
            thread = threading.Thread(target=self.on_synthesize_thread, args=(code, waveobs, regions, linelist, isotopes, abundances, atmosphere_layers, teff, logg, MH, alpha, microturbulence_vel,  macroturbulence, vsini, limb_darkening_coeff, resolution, ))
            thread.setDaemon(True)
            thread.start()

    def on_synthesize_thread(self, code, waveobs, regions, linelist, isotopes, abundances, atmosphere_layers, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, resolution):

        synth_spectrum = ispec.create_spectrum_structure(waveobs)

        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])

        error_message = None
        try:
            # waveobs is multiplied by 10.0 in order to be converted from nm to armstrongs
            synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], atmosphere_layers, teff, logg, MH, alpha, linelist, isotopes, abundances, fixed_abundances, microturbulence_vel = microturbulence_vel, macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, R=resolution, regions=regions, verbose=1, gui_queue=self.queue, code=code)
        except Exception as e:
            error_message = str(e)


        synth_spectrum.sort(order='waveobs') # Make sure it is ordered by wavelength

        # Remove atmosphere model temporary file
        self.queue.put((self.on_synthesize_finnish, [code, synth_spectrum, teff, logg, MH, alpha, microturbulence_vel, error_message], {}))

    def on_synthesize_finnish(self, code, synth_spectrum, teff, logg, MH, alpha, microturbulence_vel, error_message):
        if error_message is not None:
            msg = error_message
            title = 'Problem synthesizing spectrum'
            self.error(title, msg)
            self.operation_in_progress = False
            return
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
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        # Name: If it already exists, add a suffix
        name = self.get_name(code + "_" + str(teff) + "_" + str(logg) + "_"  + str(MH) + "_"  + str(alpha) + "_" + str(microturbulence_vel))
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

    def on_determine_abundances_from_ew(self, show_previous_results=True):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.linemasks is None:
            msg = "Lines should be fitted first."
            title = 'Lines not fitted'
            self.error(title, msg)
            self.flash_status_message("Not previous fitted lines available.")
            return

        key = "AbundancesDialog"
        if key not in self.active_spectrum.dialog:
            teff = 5771.0
            logg = 4.44
            MH = 0.00
            alpha = 0.00
            microturbulence_vel = 1.07
            self.active_spectrum.dialog[key] = AbundancesDialog(self, "Abundances determination from EW", teff, logg, MH, alpha, microturbulence_vel, self.lists, self.default_lists)
            self.active_spectrum.dialog[key].show()
        elif show_previous_results:
            self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        code = self.active_spectrum.dialog[key].results["Code"].lower()
        teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
        logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
        MH = self.active_spectrum.dialog[key].results["Metallicity [M/H]"]
        alpha = self.active_spectrum.dialog[key].results["Alpha enhancement [alpha/Fe]"]
        microturbulence_vel = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
        selected_atmosphere_models = self.active_spectrum.dialog[key].results["Model atmosphere"]
        selected_abundances = self.active_spectrum.dialog[key].results["Solar abundances"]
        abundances_file = resource_path("input/abundances/" + selected_abundances + "/stdatom.dat")
        self.active_spectrum.dialog[key].destroy()

        if teff is None or logg is None or MH is None or alpha is None or microturbulence_vel is None:
            self.flash_status_message("Bad value.")
            return


        if selected_atmosphere_models not in self.modeled_layers_pack:
            logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.modeled_layers_pack[selected_atmosphere_models] = ispec.load_modeled_layers_pack(resource_path('input/atmospheres/' + selected_atmosphere_models + '/'))

        if not ispec.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic':microturbulence_vel}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models."
            title = 'Out of the atmospheric models'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return

        # Load SPECTRUM abundances
        if not abundances_file in list(self.solar_abundances.keys()):
            self.solar_abundances[abundances_file] = ispec.read_solar_abundances(abundances_file)
        abundances = self.solar_abundances[abundances_file]

        # Prepare atmosphere model
        self.status_message("Interpolating atmosphere model...")
        atmosphere_layers = ispec.interpolate_atmosphere_layers(self.modeled_layers_pack[selected_atmosphere_models], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha})

        self.operation_in_progress = True
        self.status_message("Determining abundances...")

        thread = threading.Thread(target=self.on_determine_abundances_from_ew_thread, args=(code, atmosphere_layers, teff, logg, MH, alpha, abundances, microturbulence_vel,))
        thread.setDaemon(True)
        thread.start()


    def on_determine_abundances_from_ew_thread(self, code, atmosphere_layers, teff, logg, MH, alpha, abundances, microturbulence_vel):
        linemasks = self.active_spectrum.linemasks
        error_message = None
        try:
            spec_abund, normal_abund, x_over_h, x_over_fe = ispec.determine_abundances(atmosphere_layers, teff, logg, MH, alpha, linemasks, abundances, microturbulence_vel=microturbulence_vel, verbose=1, gui_queue=self.queue, code=code)
        except Exception as e:
            spec_abund = None
            normal_abund = None
            x_over_h = None
            x_over_fe = None
            error_message = str(e)

        self.queue.put((self.on_determine_abundances_from_ew_finnish, [spec_abund, normal_abund, x_over_h, x_over_fe, error_message], {}))

    def on_determine_abundances_from_ew_finnish(self, spec_abund, normal_abund, x_over_h, x_over_fe, error_message):
        if error_message is not None:
            msg = error_message
            title = 'Problem determining abundances'
            self.error(title, msg)
            self.operation_in_progress = False
            return

        self.flash_status_message("Abundances determined!")
        self.operation_in_progress = False

        self.active_spectrum.abundances = normal_abund

        key = "AbundancesDialog"
        #self.active_spectrum.dialog[key].register(self.active_spectrum.linemasks, spec_abund)
        self.active_spectrum.dialog[key].register(self.active_spectrum.linemasks, x_over_h, x_over_fe)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return
        else:
            # Recalculate
            self.active_spectrum.dialog[key].destroy()
            self.on_determine_abundances_from_ew(show_previous_results=False)


    def on_determine_parameters_from_ew(self, show_previous_results=True):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.linemasks is None:
            msg = "Lines should be fitted first."
            title = 'Lines not fitted'
            self.error(title, msg)
            self.flash_status_message("Not previous fitted lines available.")
            return

        key = "SolverEWDialog"
        if key not in self.active_spectrum.dialog:
            teff = 5771.0
            logg = 4.44
            MH = 0.00
            alpha = 0.00
            microturbulence_vel = 1.07
            self.active_spectrum.dialog[key] = SolverEWDialog(self, "Parameters determination from EW", teff, logg, MH, alpha, microturbulence_vel, self.lists, self.default_lists)
            self.active_spectrum.dialog[key].show()
        elif show_previous_results:
            self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        code = self.active_spectrum.dialog[key].results["Code"].lower()
        teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
        logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
        MH = self.active_spectrum.dialog[key].results["Metallicity [M/H]"]
        alpha = self.active_spectrum.dialog[key].results["Alpha enhancement [alpha/Fe]"]
        microturbulence_vel = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
        max_iterations = self.active_spectrum.dialog[key].results["Maximum number of iterations"]
        free_teff = self.active_spectrum.dialog[key].results["Free Teff"] == 1
        free_logg = self.active_spectrum.dialog[key].results["Free Log(g)"] == 1
        free_microturbulence = self.active_spectrum.dialog[key].results["Free Vmic"] == 1
        enhance_abundances = self.active_spectrum.dialog[key].results["Automatic alpha enhancement [alpha/Fe]"] == 1
        free_params = []
        if free_teff:
            free_params.append("teff")
        if free_logg:
            free_params.append("logg")
        if free_microturbulence:
            free_params.append("vmic")

        if len(free_params) == 0:
            msg = "At least one parameter should be let free"
            title = 'No free parameters'
            self.error(title, msg)
            return

        selected_atmosphere_models = self.active_spectrum.dialog[key].results["Model atmosphere"]
        selected_abundances = self.active_spectrum.dialog[key].results["Solar abundances"]
        abundances_file = resource_path("input/abundances/" + selected_abundances + "/stdatom.dat")
        self.active_spectrum.dialog[key].destroy()

        if teff is None or logg is None or MH is None or alpha is None or microturbulence_vel is None:
            self.flash_status_message("Bad value.")
            return


        if selected_atmosphere_models not in self.modeled_layers_pack:
            logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.modeled_layers_pack[selected_atmosphere_models] = ispec.load_modeled_layers_pack(resource_path('input/atmospheres/' + selected_atmosphere_models + '/'))

        if not ispec.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic':microturbulence_vel}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models."
            title = 'Out of the atmospheric models'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return

        # Load SPECTRUM abundances
        if not abundances_file in list(self.solar_abundances.keys()):
            self.solar_abundances[abundances_file] = ispec.read_solar_abundances(abundances_file)
        abundances = self.solar_abundances[abundances_file]

        self.operation_in_progress = True
        self.status_message("Determining abundances...")

        thread = threading.Thread(target=self.on_determine_parameters_from_ew_thread, args=(code, self.modeled_layers_pack[selected_atmosphere_models], abundances,  teff, logg, MH, alpha, microturbulence_vel, free_params, enhance_abundances, max_iterations))
        thread.setDaemon(True)
        thread.start()


    def on_determine_parameters_from_ew_thread(self, code, modeled_layers_pack, solar_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, free_params, enhance_abundances, max_iterations):
        linemasks = self.active_spectrum.linemasks
        # Reduced equivalent width
        # Filter too weak/strong lines
        # * Criteria presented in paper of GALA
        #efilter = np.logical_and(linemasks['ewr'] >= -5.8, linemasks['ewr'] <= -4.65)
        efilter = np.logical_and(linemasks['ewr'] >= -6.0, linemasks['ewr'] <= -4.3)
        # Filter high excitation potential lines
        # * Criteria from Eric J. Bubar "Equivalent Width Abundance Analysis In Moog"
        efilter = np.logical_and(efilter, linemasks['lower_state_eV'] <= 5.0)
        efilter = np.logical_and(efilter, linemasks['lower_state_eV'] >= 0.5)
        ## Filter also bad fits
        efilter = np.logical_and(efilter, linemasks['rms'] < 1.00)
        # no flux
        noflux = self.active_spectrum.data['flux'][linemasks['peak']] < 1.0e-10
        efilter = np.logical_and(efilter, np.logical_not(noflux))
        unfitted = linemasks['fwhm'] == 0
        efilter = np.logical_and(efilter, np.logical_not(unfitted))

        iron = np.logical_or(linemasks['element'] == "Fe 1", linemasks['element'] == "Fe 2")
        efilter = np.logical_and(efilter, iron)

        error_message = None
        try:
            results = ispec.model_spectrum_from_ew(linemasks[efilter], modeled_layers_pack, \
                                solar_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, \
                                free_params=free_params, \
                                adjust_model_metalicity=True, \
                                max_iterations=max_iterations, \
                                enhance_abundances=enhance_abundances, \
                                #outliers_detection = "robust", \
                                #outliers_weight_limit = 0.90, \
                                outliers_detection = "sigma_clipping", \
                                #sigma_level = 3, \
                                tmp_dir = None, \
                                code=code)
            params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks = results
        except Exception as e:
            params = None
            errors = None
            status = None
            x_over_h = None
            selected_x_over_h = None
            fitted_lines_params = None
            used_linemasks = None
            error_message = str(e)


        self.queue.put((self.on_determine_parameters_from_ew_finnish, [params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks, error_message], {}))

    def on_determine_parameters_from_ew_finnish(self, params, errors, status, x_over_h, selected_x_over_h, fitted_lines_params, used_linemasks, error_message):
        if error_message is not None:
            msg = error_message
            title = 'Problem determining abundances'
            self.error(title, msg)
            self.operation_in_progress = False
            return

        self.flash_status_message("Parameters determined!")
        self.operation_in_progress = False

        key = "SolverEWDialog"
        self.active_spectrum.dialog[key].register(used_linemasks, params, x_over_h, selected_x_over_h, fitted_lines_params)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            self.active_spectrum.dialog[key].destroy()
            return
        else:
            # Recalculate
            self.active_spectrum.dialog[key].destroy()
            self.on_determine_parameters_from_ew(show_previous_results=False)


    def on_determine_parameters(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        if not self.check_continuum_model_exists():
            return

        teff = 5771.0
        logg = 4.44
        MH = 0.00
        alpha = 0.00
        macroturbulence = 4.21
        vsini = 1.6
        limb_darkening_coeff = 0.6
        microturbulence_vel = 1.05
        #resolution = 47000
        resolution = 300000
        #resolution = 100000

        key = "SolverDialog"
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = SolverDialog(self, "Determine parameters", resolution, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, self.lists, self.default_lists)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        code = self.active_spectrum.dialog[key].results["Code"].lower()
        initial_teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
        initial_logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
        initial_MH = self.active_spectrum.dialog[key].results["Metallicity [M/H]"]
        initial_alpha = self.active_spectrum.dialog[key].results["Alpha enhancement [alpha/Fe]"]
        initial_vmic = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
        initial_vmac = self.active_spectrum.dialog[key].results["Macroturbulence velocity (km/s)"]
        initial_vsini = self.active_spectrum.dialog[key].results["Rotation (v sin(i)) (km/s)"]
        initial_limb_darkening_coeff = self.active_spectrum.dialog[key].results["Limb darkening coefficient"]
        initial_R = self.active_spectrum.dialog[key].results["Resolution"]
        initial_vrad = self.active_spectrum.dialog[key].results["Radial velocity"]
        selected_atmosphere_models = self.active_spectrum.dialog[key].results["Model atmosphere"]
        selected_abundances = self.active_spectrum.dialog[key].results["Solar abundances"]
        selected_linelist = self.active_spectrum.dialog[key].results["Line list"]
        element_abundance = int(self.active_spectrum.dialog[key].results["Individual abundance"].split()[0])
        #element_abundance_name = self.active_spectrum.dialog[key].results["Individual abundance"].split()[2]
        max_iterations = self.active_spectrum.dialog[key].results["Maximum number of iterations"]

        free_teff = self.active_spectrum.dialog[key].results["Free Teff"] == 1
        free_logg = self.active_spectrum.dialog[key].results["Free Log(g)"] == 1
        free_MH = self.active_spectrum.dialog[key].results["Free [M/H]"] == 1
        free_alpha = self.active_spectrum.dialog[key].results["Free [alpha/Fe]"] == 1
        free_microturbulence = self.active_spectrum.dialog[key].results["Free Vmic"] == 1
        free_macroturbulence = self.active_spectrum.dialog[key].results["Free Vmac"] == 1
        free_vsini = self.active_spectrum.dialog[key].results["Free vsin(i)"] == 1
        free_limb_darkening_coeff = self.active_spectrum.dialog[key].results["Free limb dark. coeff."] == 1
        free_resolution = self.active_spectrum.dialog[key].results["Free resolution"] == 1
        free_vrad = self.active_spectrum.dialog[key].results["Free radial velocity"] == 1
        free_element_abundance = self.active_spectrum.dialog[key].results["Free individual abundance"] == 1
        enhance_abundances = self.active_spectrum.dialog[key].results["Automatic alpha enhancement [alpha/Fe]"] == 1
        vmic_from_empirical_relation = self.active_spectrum.dialog[key].results["Automatic Vmic from empirical relation"] == 1
        vmac_from_empirical_relation = self.active_spectrum.dialog[key].results["Automatic Vmac from empirical relation"] == 1

        free_params = []
        if free_teff:
            free_params.append("teff")
        if free_logg:
            free_params.append("logg")
        if free_MH:
            free_params.append("MH")
        if free_alpha:
            free_params.append("alpha")
        if free_microturbulence:
            free_params.append("vmic")
        if free_macroturbulence:
            free_params.append("vmac")
        if free_vsini:
            free_params.append("vsini")
        if free_limb_darkening_coeff:
            free_params.append("limb_darkening_coef")
        if free_resolution:
            free_params.append("R")
        if free_vrad:
            free_params.append("vrad")
        if free_element_abundance:
            free_params.append(str(element_abundance))

        if len(free_params) == 0:
            msg = "At least one parameter should be let free"
            title = 'No free parameters'
            self.error(title, msg)
            return


        self.active_spectrum.dialog[key].destroy()

        ### Find filenames
        i = np.where(self.lists['atomic_lines']['name'] == selected_linelist)
        atomic_linelist_file = self.lists['atomic_lines']['path'][i][0]

        i = np.where(self.lists['abundances']['name'] == selected_abundances)
        abundances_file = self.lists['abundances']['path'][i][0]

        i = np.where(self.lists['atmospheres']['name'] == selected_atmosphere_models)
        atmospheres_file = os.path.dirname(self.lists['atmospheres']['path'][i][0])

        ####

        isotope_file = resource_path("input/isotopes/SPECTRUM.lst")

        if selected_atmosphere_models not in self.modeled_layers_pack:
            logging.info("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.status_message("Loading %s modeled atmospheres..." % selected_atmosphere_models)
            self.modeled_layers_pack[selected_atmosphere_models] = ispec.load_modeled_layers_pack(atmospheres_file)

        if not ispec.valid_atmosphere_target(self.modeled_layers_pack[selected_atmosphere_models], {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha, 'vmic':initial_vmic}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the atmospheric models."
            title = 'Out of the atmospheric models'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return


        # Load SPECTRUM linelist
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")
        if self.molecules is None:
            self.molecules = ispec.read_molecular_symbols(molecules_file)
        if self.chemical_elements is None:
            self.chemical_elements = ispec.read_chemical_elements(chemical_elements_file)
        if not selected_linelist in list(self.atomic_linelist.keys()):
            self.atomic_linelist[selected_linelist] = ispec.read_atomic_linelist(atomic_linelist_file)
            logging.warning("Limiting linelist to lines that have at least 0.01 depth in the Sun")
            solar = self.atomic_linelist[selected_linelist]['theoretical_depth'] >= 0.01
            self.atomic_linelist[selected_linelist] = self.atomic_linelist[selected_linelist][solar]
        linelist = self.atomic_linelist[selected_linelist]
        isotopes = ispec.read_isotope_data(isotope_file)

        # Load SPECTRUM abundances
        if not abundances_file in list(self.solar_abundances.keys()):
            self.solar_abundances[abundances_file] = ispec.read_solar_abundances(abundances_file)
        abundances = self.solar_abundances[abundances_file]

        if not free_element_abundance:
            # No fixed abundances
            free_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])
        else:
            free_abundances = np.recarray((1, ), dtype=[('code', int),('Abund', float)])
            free_abundances['code'] = element_abundance
            free_abundances['Abund'] = abundances['Abund'][abundances['code'] == int(element_abundance)] # Initial abundance


        # Consider only segments
        elements_type = "segments"
        if len(self.region_widgets[elements_type]) == 0:
            msg = "No segments present for synthetic spectrum generation."
            title = 'No segments'
            self.error(title, msg)
            return

        waveobs = self.active_spectrum.data['waveobs']
        # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectrum
        if np.min(waveobs) < 100.0 or np.max(waveobs) > 4000.0:
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
        thread = threading.Thread(target=self.on_determine_parameters_thread, args=(code, selected_atmosphere_models, linelist, isotopes, abundances, free_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, enhance_abundances, vmic_from_empirical_relation, vmac_from_empirical_relation, max_iterations))
        thread.setDaemon(True)
        thread.start()

    def on_determine_parameters_thread(self, code, selected_atmosphere_models, linelist, isotopes, abundances, free_abundances, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, enhance_abundances, vmic_from_empirical_relation, vmac_from_empirical_relation, max_iterations):
        self.__update_numpy_arrays_from_widgets("lines")
        self.__update_numpy_arrays_from_widgets("segments")

        error_message = None
        linelist_free_loggf = None
        try:
            obs_spectrum, synth_spectrum, params, errors, free_abundances, loggf_found, status, stats_linemasks = ispec.model_spectrum(self.active_spectrum.data, self.active_spectrum.continuum_model, self.modeled_layers_pack[selected_atmosphere_models], linelist, isotopes, abundances, free_abundances, linelist_free_loggf, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=self.regions['segments'], linemasks=self.regions['lines'], max_iterations=max_iterations, code=code, \
            enhance_abundances=enhance_abundances, \
            use_errors = True, \
            vmic_from_empirical_relation = vmic_from_empirical_relation, \
            vmac_from_empirical_relation = vmac_from_empirical_relation \
            )
        except Exception as e:
            raise
            obs_spectrum = None
            synth_spectrum = None
            params = None
            errors = None
            status = None
            stats_linemasks = None
            error_message = str(e)
        self.queue.put((self.on_determine_parameters_finnish, [obs_spectrum, synth_spectrum, params, errors, status, stats_linemasks, error_message], {}))

    def on_determine_parameters_finnish(self, obs_spectrum, synth_spectrum, params, errors, status, stats_linemasks, error_message):
        if error_message is not None:
            msg = error_message
            title = 'Problem synthesizing spectrum'
            self.error(title, msg)
            self.operation_in_progress = False
            return
        # Name: If it already exists, add a suffix
        base_name = "%.1f_%.2f_%.2f_%.2f_%.1f_%.1f_%.1f_%.1f_%.0f" % (params['teff'], params['logg'], params['MH'], params['alpha'], params['vmic'], params['vmac'], params['vsini'], params['limb_darkening_coeff'], params['R'])

        analysed_spectrum = self.active_spectrum

        ########## Synthetic
        name = self.get_name(base_name)
        color = self.get_color()
        ### Add a new spectrum but do not make it active to avoid confusions
        self.active_spectrum = Spectrum(synth_spectrum, name, color=color)
        self.spectra.append(self.active_spectrum)
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        ########## Observed
        #name = self.get_name(base_name + "_obs")
        #color = self.get_color()
        #### Add a new spectrum but do not make it active to avoid confusions
        #self.active_spectrum = Spectrum(obs_spectrum, name, color=color)
        #self.spectra.append(self.active_spectrum)
        #self.active_spectrum.not_saved = True
        #self.draw_active_spectrum()

        # Remove "[A]  " from spectrum name (legend) if it exists
        if self.active_spectrum is not None and self.active_spectrum.plot_id is not None:
            self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        self.active_spectrum = analysed_spectrum
        self.draw_active_spectrum()

        self.update_title()
        self.update_menu_active_spectrum()
        self.update_scale()

        self.canvas.draw()

        self.operation_in_progress = False
        self.flash_status_message("Parameters determined!")

    def on_determine_parameters_with_grid(self):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return
        if not self.check_continuum_model_exists():
            return

        teff = 5771.0
        logg = 4.44
        MH = 0.00
        alpha = 0.00
        macroturbulence = 4.21
        vsini = 1.6
        limb_darkening_coeff = 0.6
        microturbulence_vel = 1.05
        #resolution = 47000
        resolution = 300000
        #resolution = 100000

        key = "InterpolateSolverDialog"
        if key not in self.active_spectrum.dialog:
            self.active_spectrum.dialog[key] = InterpolateSolverDialog(self, "Determine parameters with grid", resolution, teff, logg, MH, alpha, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, self.lists, self.default_lists)
        self.active_spectrum.dialog[key].show()

        if self.active_spectrum.dialog[key].results is None:
            # Cancel
            self.active_spectrum.dialog[key].destroy()
            return

        selected_grid = self.active_spectrum.dialog[key].results["Grid"]
        initial_teff = self.active_spectrum.dialog[key].results["Effective temperature (K)"]
        initial_logg = self.active_spectrum.dialog[key].results["Surface gravity (log g)"]
        initial_MH = self.active_spectrum.dialog[key].results["Metallicity [M/H]"]
        initial_alpha = self.active_spectrum.dialog[key].results["Alpha enhancement [alpha/Fe]"]
        initial_vmic = self.active_spectrum.dialog[key].results["Microturbulence velocity (km/s)"]
        initial_vmac = self.active_spectrum.dialog[key].results["Macroturbulence velocity (km/s)"]
        initial_vsini = self.active_spectrum.dialog[key].results["Rotation (v sin(i)) (km/s)"]
        initial_limb_darkening_coeff = self.active_spectrum.dialog[key].results["Limb darkening coefficient"]
        initial_R = self.active_spectrum.dialog[key].results["Resolution"]
        initial_vrad = self.active_spectrum.dialog[key].results["Radial velocity"]
        max_iterations = self.active_spectrum.dialog[key].results["Maximum number of iterations"]

        free_teff = self.active_spectrum.dialog[key].results["Free Teff"] == 1
        free_logg = self.active_spectrum.dialog[key].results["Free Log(g)"] == 1
        free_MH = self.active_spectrum.dialog[key].results["Free [M/H]"] == 1
        free_alpha = self.active_spectrum.dialog[key].results["Free [alpha/Fe]"] == 1
        free_microturbulence = self.active_spectrum.dialog[key].results["Free Vmic"] == 1
        free_macroturbulence = self.active_spectrum.dialog[key].results["Free Vmac"] == 1
        free_vsini = self.active_spectrum.dialog[key].results["Free vsin(i)"] == 1
        free_limb_darkening_coeff = self.active_spectrum.dialog[key].results["Free limb dark. coeff."] == 1
        free_resolution = self.active_spectrum.dialog[key].results["Free resolution"] == 1
        free_vrad = self.active_spectrum.dialog[key].results["Free radial velocity"] == 1
        vmac_from_empirical_relation = self.active_spectrum.dialog[key].results["Automatic Vmac from empirical relation"] == 1

        free_params = []
        if free_teff:
            free_params.append("teff")
        if free_logg:
            free_params.append("logg")
        if free_MH:
            free_params.append("MH")
        if free_alpha:
            free_params.append("alpha")
        if free_microturbulence:
            free_params.append("vmic")
        if free_macroturbulence:
            free_params.append("vmac")
        if free_vsini:
            free_params.append("vsini")
        if free_limb_darkening_coeff:
            free_params.append("limb_darkening_coef")
        if free_resolution:
            free_params.append("R")
        if free_vrad:
            free_params.append("vrad")

        if len(free_params) == 0:
            msg = "At least one parameter should be let free"
            title = 'No free parameters'
            self.error(title, msg)
            return

        self.active_spectrum.dialog[key].destroy()

        ### Find filenames
        i = np.where(self.lists['grid']['name'] == selected_grid)
        grid_dirname = os.path.dirname(self.lists['grid']['path'][i][0])

        ####

        if selected_grid not in self.grid:
            logging.info("Loading grid %s..." % selected_grid)
            self.status_message("Loading grid %s..." % selected_grid)
            self.grid[selected_grid] = ispec.load_spectral_grid(grid_dirname)

        if not ispec.valid_interpolated_spectrum_target(self.grid[selected_grid], {'teff':initial_teff, 'logg':initial_logg, 'MH':initial_MH, 'alpha':initial_alpha, 'vmic':initial_vmic}):
            msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of the grid limits."
            title = 'Out of the grid limits'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return

        # Consider only segments
        elements_type = "segments"
        if len(self.region_widgets[elements_type]) == 0:
            msg = "No segments present for synthetic spectrum generation."
            title = 'No segments'
            self.error(title, msg)
            return

        waveobs = self.active_spectrum.data['waveobs']
        total_points = len(waveobs)

        if total_points < 2:
            msg = "Wavelength range too narrow."
            title = 'Wavelength range'
            self.error(title, msg)
            self.flash_status_message("Bad values.")
            return

        self.operation_in_progress = True
        self.status_message("Determining parameters with grid...")
        self.update_progress(10)
        code = "grid"
        thread = threading.Thread(target=self.on_determine_parameters_with_grid_thread, args=(code, selected_grid, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, vmac_from_empirical_relation, max_iterations))
        thread.setDaemon(True)
        thread.start()

    def on_determine_parameters_with_grid_thread(self, code, selected_grid, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, vmac_from_empirical_relation, max_iterations):
        self.__update_numpy_arrays_from_widgets("lines")
        self.__update_numpy_arrays_from_widgets("segments")

        # No fixed abundances
        fixed_abundances = np.recarray((0, ), dtype=[('code', int),('Abund', float)])
        modeled_layers_pack = None
        linelist = None
        isotopes = None
        abundances = None
        error_message = None
        linelist_free_loggf = None
        free_abundances = None
        try:
            obs_spectrum, synth_spectrum, params, errors, free_abundances, loggf_found, status, stats_linemasks = ispec.model_spectrum(self.active_spectrum.data, self.active_spectrum.continuum_model, modeled_layers_pack, linelist, isotopes, abundances, free_abundances, linelist_free_loggf, initial_teff, initial_logg, initial_MH, initial_alpha, initial_vmic, initial_vmac, initial_vsini, initial_limb_darkening_coeff, initial_R, initial_vrad, free_params, segments=self.regions['segments'], linemasks=self.regions['lines'], max_iterations=max_iterations, code=code, \
            vmic_from_empirical_relation=False, \
            vmac_from_empirical_relation=vmac_from_empirical_relation, \
            grid=self.grid[selected_grid] \
            )
        except Exception as e:
            #raise
            obs_spectrum = None
            synth_spectrum = None
            params = None
            errors = None
            status = None
            stats_linemasks = None
            error_message = str(e)
        self.queue.put((self.on_determine_parameters_finnish, [obs_spectrum, synth_spectrum, params, errors, status, stats_linemasks, error_message], {}))

