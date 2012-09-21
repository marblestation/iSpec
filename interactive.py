#!/usr/bin/env python
#
#    This file is part of Spectra Visual Editor (SVE).
#    Copyright 2011-2012 Sergi Blanco Cuaresma - http://www.marblestation.com
#
#    SVE is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SVE is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with SVE. If not, see <http://www.gnu.org/licenses/>.
#
#################
# Run with ipython -pdb -c "%run interactive.py"
#################
#import ipdb

import os
import sys


import pprint
import random
import wx
import numpy as np
import numpy.lib.recfunctions as rfn # Extra functions

import threading
import getopt

import logging

# The recommended way to use wx with mpl is with the WXAgg backend.
import matplotlib
matplotlib.use('WXAgg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar
from matplotlib.ticker import ScalarFormatter

import sve
from dialogs import *

## PyInstaller resource access
def resource_path(relative):
    if getattr(sys, 'frozen', None):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(__file__)
    return os.path.join(basedir, relative)

class CustomizableRegion:
    min_width = 0.002 # nm
    mark = None

    def __init__(self, frame, element_type, axvspan, mark=None, note=None):
        self.frame = frame
        self.axvspan = axvspan
        self.element_type = element_type
        self.press = None
        self.original_facecolor = self.axvspan.get_facecolor()
        self.original_edgecolor = self.axvspan.get_edgecolor()
        ## Line region specific properties:
        self.mark = mark
        self.note = note
        if mark != None:
            self.mark_position = self.mark.get_xdata()
        else:
            self.mark_position = None

        # Fit line properties, dictionary of spectrum (to vinculate to different spectrum):
        self.line_plot_id = {}
        self.line_model = {}
        self.line_extra = {}

    def connect(self):
        # Connect to all the events
        self.cid_press = self.axvspan.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cid_release = self.axvspan.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cid_motion = self.axvspan.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)


    def update_size(self, event):
        # Update the position of the edge left or right
        xy = self.axvspan.get_xy()

        # If left button clicked, modify left edge
        if event.button == 1:
            # Check condition to allow modification
            if self.mark != None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata < self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True

            # Do not modify if the region will become slimmer than...
            big_enough = xy[2,0] - event.xdata > self.min_width

            if big_enough and compatible_with_mark:
                xy[0,0] = event.xdata
                xy[1,0] = event.xdata
                xy[4,0] = event.xdata
                self.frame.status_message("Moving left edge to %.4f" % xy[0,0])
                self.frame.canvas.draw()
        elif event.button == 3:
            # Check condition to allow modification
            if self.mark != None:
                # Mark should be inside the region
                compatible_with_mark = event.xdata > self.mark.get_xdata()[0]
            else:
                # There is no mark, so any position is good
                compatible_with_mark = True

            # Do not modify if the region will become slimmer than...
            big_enough = event.xdata - xy[0,0] > self.min_width

            if big_enough and compatible_with_mark:
                xy[2,0] = event.xdata
                xy[3,0] = event.xdata
                self.frame.status_message("Moving right edge to %.4f" % xy[2,0])
                self.frame.canvas.draw()

    def update_mark(self, event):
        # Position of the edge left or right
        xy = self.axvspan.get_xy()

        # If left button clicked, modify mark
        if event.button == 1 and self.mark != None:
            x = self.mark.get_xdata()

            inside_region = (event.xdata > xy[0,0]) and (event.xdata < xy[2,0])
            if inside_region:
                x[0] = event.xdata
                x[1] = event.xdata
                self.mark.set_xdata(x)
                self.mark_position = self.mark.get_xdata()
                if self.note != None:
                    self.note.xy = (event.xdata, 1)
                self.frame.status_message("Moving mark to %.4f" % x[0])
                self.frame.canvas.draw()

    def update_mark_note(self, event):
        if self.note != None:
            note_text = self.note.get_text()
        else:
            note_text = ""

        note_text = self.frame.ask_value('Note for the new line region:', 'Note', note_text)
        self.update_mark_note_text(note_text)
        self.frame.canvas.draw()

    def update_mark_note_text(self, note_text):
        if note_text != None and note_text != "":
            # Set
            if self.note == None:
                # Position of the mark
                x = self.mark_position[0]

                # New
                self.note = self.frame.axes.annotate(note_text, xy=(x, 1),  xycoords=("data", 'axes fraction'),
                    xytext=(-10, 20), textcoords='offset points',
                    size=8,
                    bbox=dict(boxstyle="round", fc="0.8"),
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                    edgecolor='black'),
                    horizontalalignment='right', verticalalignment='top',
                    annotation_clip=True,
                    )
            else:
                # Update
                self.note.set_text(note_text)
        elif note_text == "" and self.note != None:
            # Remove
            self.note.set_visible(False)
            self.note = None

    def show_mark(self):
        if self.mark == None:
            self.mark = self.frame.axes.axvline(x = self.mark_position[0], linewidth=1, color='orange')
            pass

    def hide_mark(self):
        if self.mark != None:
            self.mark.remove()
            self.mark = None


    def on_press(self, event):
        # Validate that the click is on this axis/object
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return
        contains, attrd = self.axvspan.contains(event)
        if not contains: return
        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return

        if self.frame.operation_in_progress:
            return

        # When regions overlap two or more can receive the click event, so
        # let's use a lock to allow the modification of one of them
        if self.frame.lock.acquire(False):
            if self.frame.action == "Remove":
                # The real removal is in the "on_release" event
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    self.press = event.button, event.x, event.xdata
                else:
                    self.press = event.button, event.x, event.xdata
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Modify":
                ## Modify region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata

                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    if event.button == 1:
                        # Left button, modify mark position
                        self.update_mark(event)
                    elif event.button == 2:
                        pass
                        # Modification of the note is managed on_release
                        # because we cannot show modal windows asking information
                        # while on_press or the mouse stops working
                else:
                    self.update_size(event)
                self.frame.regions_changed(self.element_type)
                # Do not free the lock until the user releases the click
            elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                ## Create note but handel it on_release
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata
            elif self.frame.action == "Stats":
                ## Statistics about the region
                self.axvspan.set_color('red')
                self.axvspan.set_alpha(0.5)
                self.frame.canvas.draw()
                self.press = event.button, event.x, event.xdata


    def disconnect_and_remove(self):
        #self.hide()
        #self.frame.canvas.draw()
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_press)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_release)
        self.axvspan.figure.canvas.mpl_disconnect(self.cid_motion)
        self.axvspan.remove()
        if self.mark != None:
            self.mark.remove()
        if self.note != None:
            self.note.remove()

    def hide(self):
        self.axvspan.set_visible(False)
        if self.mark != None:
            self.mark.set_visible(False)
        if self.note != None:
            self.note.set_visible(False)
        if self.line_plot_id != None:
            for plot_id in self.line_plot_id.values():
                if plot_id != None:
                    self.frame.axes.lines.remove(plot_id)


    def on_release(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return

        # Validate it is not in PAN or ZOOM mode
        if event.inaxes != None and event.inaxes.get_navigate_mode() != None: return

        # This element is the kind of element that is selected to be modified?
        if self.frame.elements != self.element_type: return
        # If the action is "create", this should be managed by the frame and not individual elements
        if self.frame.action == "Create" and not (self.frame.elements == "lines" and self.frame.subelements == "marks"): return

        if self.press != None:
            # In modification mode, if it is the current selected widget
            self.press = None
            if self.frame.action == "Remove":
                self.frame.lock.release() # Release now because the next commands will
                if self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Remove note
                    if self.note != None:
                        self.note.set_visible(False)
                        self.note.remove()
                        self.note = None
                        self.frame.canvas.draw()
                else:
                    # Remove from the figure now (except for line marks)
                    #  (if we do it on_press, this event would not be trigered and the lock not released)
                    self.frame.region_widgets[self.element_type].remove(self)
                    self.frame.regions_changed(self.element_type)
                    self.frame.flash_status_message("Removed region from " + "%.4f" % self.axvspan.get_xy()[0,0] + " to " + "%.4f" % self.axvspan.get_xy()[2,0])
                    self.disconnect_and_remove()
                    self.hide()
                    self.frame.canvas.draw()
            else:
                if self.frame.action == "Modify":
                    # If it was a right click when elements is lines and subelements marks
                    # => Change note
                    if event.button == 3 and self.frame.elements == "lines" and self.frame.subelements == "marks":
                        # Right button, modify note
                        if self.frame.action == "Modify" or self.frame.action == "Create":
                            self.update_mark_note(event)
                        else:
                            # Note removal is managed on_press
                            pass
                elif self.frame.action == "Create" and self.frame.elements == "lines" and self.frame.subelements == "marks":
                    # Create a note (left or right click)
                    if (event.button == 1 or event.button == 3) and self.note == None:
                        self.update_mark_note(event)
                elif self.frame.action == "Stats" and self.frame.lock.locked():
                    self.frame.update_stats(self)

                self.axvspan.set_facecolor(self.original_facecolor)
                self.axvspan.set_edgecolor(self.original_edgecolor)
                self.axvspan.set_alpha(0.30)
                self.frame.status_message("")
                self.frame.lock.release()
                self.frame.canvas.draw()

    def on_motion(self, event):
        # Validate that the object has been pressed and the click is on this axis
        if self.press is None: return
        if event.inaxes != self.axvspan.axes: return
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes.get_navigate_mode() != None: return
        if self.frame.action == "Stats": return

#        button, x, xpress = self.press
        if self.frame.action == "Modify":
            if self.frame.subelements == "marks":
                self.update_mark(event)
            else:
                self.update_size(event)


    def get_note_text(self):
        note_text = ""
        if self.element_type == "lines" and self.note != None:
            note_text = self.note.get_text()
        return note_text

    def get_wave_peak(self):
        if self.element_type == "lines":
            return self.mark.get_xdata()[0]
        else:
            return None

    def get_wave_base(self):
        return self.axvspan.get_xy()[0,0]

    def get_wave_top(self):
        return self.axvspan.get_xy()[2,0]


class Spectrum():
    def __init__(self, data, name, color='b', not_saved=False, plot_id=None, continuum_model=None, continuum_data=None, continuum_plot_id=None, path=None):
        self.data = data
        self.name = name
        self.color = color
        self.not_saved = not_saved
        self.plot_id = plot_id
        self.continuum_model = continuum_model
        self.continuum_data = continuum_data
        self.continuum_plot_id = continuum_plot_id
        self.errors_plot_id1 = None
        self.errors_plot_id2 = None
        self.velocity_atomic = 0.0 # Radial velocity (km/s)
        self.velocity_telluric = 0.0 # Barycentric velocity (km/s)
        self.path = path
        # Resolution power from the velocity profile relative to atomic and telluric lines
        self.resolution_atomic = 0.0
        self.resolution_telluric = 0.0
        self.velocity_profile_atomic_xcoord = None
        self.velocity_profile_atomic_fluxes = None
        self.velocity_profile_atomic_models = None
        self.velocity_profile_atomic_num_used_lines = None
        self.velocity_profile_atomic_rv_step = None
        self.velocity_profile_atomic_snr = 0.0
        self.velocity_profile_telluric_xcoord = None
        self.velocity_profile_telluric_fluxes = None
        self.velocity_profile_telluric_models = None
        self.velocity_profile_telluric_num_used_lines = None
        self.velocity_profile_telluric_rv_step = None
        self.velocity_profile_telluric_fwhm_correction = 0.0
        self.snr = None

class SVETaskBarIcon(wx.TaskBarIcon):
    TBMENU_RESTORE = wx.NewId()
    TBMENU_CLOSE   = wx.NewId()
    TBMENU_CHANGE  = wx.NewId()
    TBMENU_REMOVE  = wx.NewId()

    def __init__(self, frame):
        #version = wx.version().split()[0].split(".")
        #version = map(int, version)
        #if version[0] == 2 and version[1] == 9 and version[2] > 2:
            #wx.TaskBarIcon.__init__(self, iconType=wx.TBI_DOCK)
        #else:
            #wx.TaskBarIcon.__init__(self)
        try:
            # Only from wxPython 2.9.2.0, the argument iconType exists
            wx.TaskBarIcon.__init__(self, iconType=wx.TBI_DOCK)
        except Exception:
            wx.TaskBarIcon.__init__(self)

        self.frame = frame

        # Set the image
        icon = wx.Icon(resource_path("images/SVE.ico"), wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon, "SVE")
        self.imgidx = 1

        # bind some events
        self.Bind(wx.EVT_TASKBAR_LEFT_DCLICK, self.OnTaskBarActivate)
        self.Bind(wx.EVT_MENU, self.OnTaskBarActivate, id=self.TBMENU_RESTORE)
        self.Bind(wx.EVT_MENU, self.OnTaskBarClose, id=self.TBMENU_CLOSE)

    def CreatePopupMenu(self):
        """
        This method is called by the base class when it needs to popup
        the menu for the default EVT_RIGHT_DOWN event.  Just create
        the menu how you want it and return it from this function,
        the base class takes care of the rest.
        """
        menu = wx.Menu()
        #menu.Append(self.TBMENU_RESTORE, "Restore")
        #menu.Append(self.TBMENU_CLOSE,   "Close")
        return menu


    def MakeIcon(self, img):
        """
        The various platforms have different requirements for the
        icon size...
        """
        if "wxMSW" in wx.PlatformInfo:
            img = img.Scale(16, 16)
        elif "wxGTK" in wx.PlatformInfo:
            img = img.Scale(22, 22)
        # wxMac can be any size upto 128x128, so leave the source img alone....
        icon = wx.IconFromBitmap(img.ConvertToBitmap() )
        return icon


    def OnTaskBarActivate(self, evt):
        if self.frame.IsIconized():
            self.frame.Iconize(False)
        if not self.frame.IsShown():
            self.frame.Show(True)
        self.frame.Raise()


    def OnTaskBarClose(self, evt):
        wx.CallAfter(self.frame.Close)


class SpectraFrame(wx.Frame):
    title = 'Spectra Visual Editor '

    def ipython(self):
        import IPython
        self.embedshell = IPython.Shell.IPShellEmbed(argv=['--pdb'])
        self.embedshell()

    # regions should be a dictionary with 'continuum', 'lines' and 'segments' keys
    # filenames should be a dictionary with 'spectra', 'continuum', 'lines' and 'segments' keys
    def __init__(self, spectra=None, regions=None, filenames=None):
        self.embedshell = None
        self.ipython_thread = None
        #self.ipython_thread = threading.Thread(target=self.ipython)
        #self.ipython_thread.setDaemon(True)
        #self.ipython_thread.start()

        self.spectra_colors = ('#0000FF', '#A52A2A', '#A020F0', '#34764A', '#000000', '#90EE90', '#FFA500', '#1E90FF',   '#FFC0CB', '#7F7F7F', '#00FF00',)

        self.spectra = []
        self.active_spectrum = None
        for i in np.arange(len(spectra)):
            path = filenames["spectra"][i]
            name = path.split('/')[-1]
            name = self.get_name(name) # If it already exists, add a suffix
            color = self.get_color()
            self.active_spectrum = Spectrum(spectra[i], name, path=path, color=color)
            self.spectra.append(self.active_spectrum)


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

        self.region_widgets = {} # It is important not to lose the reference to the regions or the garbage
        self.region_widgets['continuum'] = []
        self.region_widgets['lines'] = []
        self.region_widgets['segments'] = []

        self.action = "Stats"
        self.elements = "continuum"
        self.subelements = None
        self.lock = threading.Lock()
        self.not_saved = {}
        self.not_saved["continuum"] = False
        self.not_saved["lines"] = False
        self.not_saved["segments"] = False
        self.velocity_telluric_lower_limit = -100 # km/s
        self.velocity_telluric_upper_limit = 100 # km/s
        self.velocity_telluric_step = 0.5 # km/s
        self.velocity_atomic_lower_limit = -200 # km/s
        self.velocity_atomic_upper_limit = 200 # km/s
        self.velocity_atomic_step = 1.0 # km/s
        self.linelist_atomic = None
        self.linelist_telluric = None
        self.modeled_layers_pack = None # Synthesize spectrum (atmospheric models)
        self.find_continuum_regions_wave_step = 0.05
        self.find_continuum_regions_sigma = 0.001
        self.find_continuum_regions_max_continuum_diff = 1.0
        self.find_lines_min_depth = 0.05 # (% of the continuum)
        self.find_lines_max_depth = 1.00 # (% of the continuum)
        self.show_errors = False

        # Barycentric velocity determination (default params)
        self.day = 15
        self.month = 2
        self.year = 2012
        self.hours = 0
        self.minutes = 0
        self.seconds = 0
        self.ra_hours = 19
        self.ra_minutes = 50
        self.ra_seconds = 46.99
        self.dec_degrees = 8
        self.dec_minutes = 52
        self.dec_seconds = 5.96
        self.barycentric_vel = 0.0 # km/s

        self.operation_in_progress = False

        wx.Frame.__init__(self, None, -1, self.title)
        self.tbicon = SVETaskBarIcon(self)
        self.icon = wx.Icon(resource_path("images/SVE.ico"), wx.BITMAP_TYPE_ICO)
        wx.Frame.SetIcon(self, self.icon)

        self.Bind(wx.EVT_CLOSE, self.on_close)
        self.create_menu()
        self.create_status_bar()
        self.create_main_panel()
        self.timeroff = wx.Timer(self)
        self.draw_figure_for_first_time()

    ####################################################################
    def error(self, title, msg):
        dlg_error = wx.MessageDialog(self, msg, title, wx.OK|wx.ICON_ERROR)
        dlg_error.ShowModal()
        dlg_error.Destroy()

    def info(self, title, msg):
        dlg_error = wx.MessageDialog(self, msg, title, wx.OK|wx.ICON_INFORMATION)
        dlg_error.ShowModal()
        dlg_error.Destroy()

    def question(self, title, msg):
        dlg_question = wx.MessageDialog(self, msg, title, wx.YES|wx.NO|wx.ICON_QUESTION)
        answer_yes = (dlg_question.ShowModal() == wx.ID_YES)
        dlg_question.Destroy()
        if not answer_yes:
            self.flash_status_message("Discarded")
        return answer_yes

    def ask_value(self, text, title, default_value):
        dlg_question = wx.TextEntryDialog(self, text, title)

        dlg_question.SetValue(default_value)
        if dlg_question.ShowModal() == wx.ID_OK:
            response = dlg_question.GetValue()
        else:
            response = None

        dlg_question.Destroy()
        return response
    ####################################################################

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
    ####################################################################

    # Check if exists a spectrum with that name, in that case add a suffix
    def get_name(self, name):
        num_repeated = 0
        max_num = 0
        for spec in self.spectra:
            if spec.name.startswith(name):
                try:
                    # Does it has already a suffix?
                    num = int(spec.name.split("-")[-1])
                    if num > max_num:
                        max_num = num
                except ValueError as e:
                    pass
                num_repeated += 1

        if num_repeated > 0: # There are repeated names
            name = name + "-" + str(max_num+1) # Add identificator number + 1
        return name


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
    ####################################################################

    def create_menu(self):
        self.menubar = wx.MenuBar()

        self.spectrum_function_items = []

        menu_file = wx.Menu()
        m_expt = menu_file.Append(-1, "&Open spectra\tCtrl-O", "Open spectra file")
        self.Bind(wx.EVT_MENU, self.on_open_spectra, m_expt)
        m_expt = menu_file.Append(-1, "Open co&ntinuum regions\tCtrl-N", "Open continuum regions file")
        self.Bind(wx.EVT_MENU, self.on_open_continuum, m_expt)
        m_expt = menu_file.Append(-1, "Open l&ine regions\tCtrl-I", "Open line regions file")
        self.Bind(wx.EVT_MENU, self.on_open_lines, m_expt)
        m_expt = menu_file.Append(-1, "Open se&gments\tCtrl-G", "Open segments file")
        self.Bind(wx.EVT_MENU, self.on_open_segments, m_expt)
        menu_file.AppendSeparator()
        m_expt = menu_file.Append(-1, "&Save plot image as\tCtrl-S", "Save plot image to file")
        self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
        m_expt = menu_file.Append(-1, "Save s&pectrum as\tCtrl-P", "Save active spectrum to file")
        self.Bind(wx.EVT_MENU, self.on_save_spectrum, m_expt)
        self.spectrum_function_items.append(m_expt)
        m_expt = menu_file.Append(-1, "Save &continuum regions as\tCtrl-C", "Save continuum regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_continuum_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save &line regions as\tCtrl-L", "Save line regions to file")
        self.Bind(wx.EVT_MENU, self.on_save_line_regions, m_expt)
        m_expt = menu_file.Append(-1, "Save s&egments as\tCtrl-E", "Save segments to file")
        self.Bind(wx.EVT_MENU, self.on_save_segments, m_expt)
        menu_file.AppendSeparator()
        m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
        self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

        menu_edit = wx.Menu()
        self.menu_active_spectrum = wx.Menu()
        m_expt = menu_edit.AppendMenu(-1, 'Select spectrum', self.menu_active_spectrum)
        self.spectrum_function_items.append(m_expt)

        m_close_spectrum = menu_edit.Append(-1, "Close spectrum", "Close current active spectrum")
        self.Bind(wx.EVT_MENU, self.on_close_spectrum, m_close_spectrum)
        self.spectrum_function_items.append(m_close_spectrum)
        m_close_spectra = menu_edit.Append(-1, "Close all spectra", "Close all loaded spectra")
        self.Bind(wx.EVT_MENU, self.on_close_all_spectra, m_close_spectra)
        self.spectrum_function_items.append(m_close_spectra)
        menu_edit.AppendSeparator()

        m_fit_continuum = menu_edit.Append(-1, "Fit &continuum", "Fit continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_continuum, m_fit_continuum)
        self.spectrum_function_items.append(m_fit_continuum)
        m_fit_lines = menu_edit.Append(-1, "Fit &lines", "Fit lines using the fitted continuum")
        self.Bind(wx.EVT_MENU, self.on_fit_lines, m_fit_lines)
        self.spectrum_function_items.append(m_fit_lines)
        menu_edit.AppendSeparator()

        #####
        menu_clear = wx.Menu()
        m_expt = menu_edit.AppendMenu(-1, 'Clear...', menu_clear)

        m_remove_fitted_continuum = menu_clear.Append(-1, "Fitted continuum", "Remove the fitted continuum")
        self.Bind(wx.EVT_MENU, self.on_remove_fitted_continuum, m_remove_fitted_continuum)
        self.spectrum_function_items.append(m_remove_fitted_continuum)

        m_remove_fitted_lines = menu_clear.Append(-1, "Fitted lines", "Remove fitted lines")
        self.Bind(wx.EVT_MENU, self.on_remove_fitted_lines, m_remove_fitted_lines)
        self.spectrum_function_items.append(m_remove_fitted_lines)

        m_remove_continuum_regions = menu_clear.Append(-1, "Continuum regions", "Clear continuum regions")
        self.Bind(wx.EVT_MENU, self.on_remove_continuum_regions, m_remove_continuum_regions)

        m_remove_line_masks = menu_clear.Append(-1, "Line masks", "Clear line masks")
        self.Bind(wx.EVT_MENU, self.on_remove_line_masks, m_remove_line_masks)

        m_remove_segments = menu_clear.Append(-1, "Segments", "Clear segments")
        self.Bind(wx.EVT_MENU, self.on_remove_segments, m_remove_segments)
        menu_edit.AppendSeparator()

        m_find_continuum = menu_edit.Append(-1, "&Find continuum regions", "Find continuum regions")
        self.Bind(wx.EVT_MENU, self.on_find_continuum, m_find_continuum)
        self.spectrum_function_items.append(m_find_continuum)

        m_find_lines = menu_edit.Append(-1, "&Find line masks", "Find line masks")
        self.Bind(wx.EVT_MENU, self.on_find_lines, m_find_lines)
        self.spectrum_function_items.append(m_find_lines)
        menu_edit.AppendSeparator()

        #####
        menu_determine_velocity = wx.Menu()
        m_expt = menu_edit.AppendMenu(-1, 'Determine velocity relative to...', menu_determine_velocity)
        self.spectrum_function_items.append(m_expt)

        m_determine_rv = menu_determine_velocity.Append(-1, "&Atomic lines (radial velocity)", "Determine radial velocity using atomic lines")
        self.Bind(wx.EVT_MENU, self.on_determine_velocity_atomic, m_determine_rv)

        m_determine_rv = menu_determine_velocity.Append(-1, "&Telluric lines  (barycentric velocity)", "Determine radial/barycentric velocity using telluric lines")
        self.Bind(wx.EVT_MENU, self.on_determine_velocity_telluric, m_determine_rv)
        #####
        #####
        menu_correct_velocity = wx.Menu()
        m_expt = menu_edit.AppendMenu(-1, 'Correct velocity...', menu_correct_velocity)
        self.spectrum_function_items.append(m_expt)

        m_correct_rv = menu_correct_velocity.Append(-1, "&Radial velocity", "Correct spectrum by using its radial velocity")
        self.Bind(wx.EVT_MENU, self.on_correct_rv, m_correct_rv)
        self.spectrum_function_items.append(m_correct_rv)

        m_correct_barycentric_vel = menu_correct_velocity.Append(-1, "Barycentric velocity", "Correct spectrum by using its radial velocity")
        self.Bind(wx.EVT_MENU, self.on_correct_barycentric_vel, m_correct_barycentric_vel)
        self.spectrum_function_items.append(m_correct_barycentric_vel)
        #####

        m_determine_barycentric_vel = menu_edit.Append(-1, "&Calculate barycentric velocity", "Calculate baricentryc velocity")
        self.Bind(wx.EVT_MENU, self.on_determine_barycentric_vel, m_determine_barycentric_vel)


        menu_edit.AppendSeparator()

        m_estimate_snr = menu_edit.Append(-1, "&Estimate global SNR", "Estimate Signal-to-Noise Ratio using the whole spectrum")
        self.Bind(wx.EVT_MENU, self.on_estimate_snr, m_estimate_snr)
        self.spectrum_function_items.append(m_estimate_snr)

        m_estimate_errors = menu_edit.Append(-1, "Estimate errors based on SNR", "Estimate errors by using the Signal-to-Noise Ratio")
        self.Bind(wx.EVT_MENU, self.on_estimate_errors, m_estimate_errors)
        self.spectrum_function_items.append(m_estimate_errors)

        self.m_show_errors = menu_edit.Append(-1, "Show errors in plot", "Show errors in plot", kind=wx.ITEM_CHECK)
        self.Bind(wx.EVT_MENU, self.on_show_errors, self.m_show_errors)
        self.spectrum_function_items.append(self.m_show_errors)
        #self.m_show_errors.Check()
        #self.m_show_errors.Check(False)

        menu_edit.AppendSeparator()

        m_degrade_resolution = menu_edit.Append(-1, "Degrade resolution", "Degread spectrum resolution")
        self.Bind(wx.EVT_MENU, self.on_degrade_resolution, m_degrade_resolution)
        self.spectrum_function_items.append(m_degrade_resolution)
        m_normalize_spectrum = menu_edit.Append(-1, "Continuum normalization", "Normalize spectrum")
        self.Bind(wx.EVT_MENU, self.on_continuum_normalization, m_normalize_spectrum)
        self.spectrum_function_items.append(m_normalize_spectrum)
        m_clean_spectrum = menu_edit.Append(-1, "Clean fluxes and errors", "Filter out measurements with fluxes and errors out of range")
        self.Bind(wx.EVT_MENU, self.on_clean_spectrum, m_clean_spectrum)
        self.spectrum_function_items.append(m_clean_spectrum)
        m_cut_spectrum = menu_edit.Append(-1, "Wavelength range reduction", "Reduce the wavelength range")
        self.Bind(wx.EVT_MENU, self.on_cut_spectrum, m_cut_spectrum)
        self.spectrum_function_items.append(m_cut_spectrum)
        m_convert_to_nm = menu_edit.Append(-1, "Convert to nanometers", "Divide wavelength by 10")
        self.Bind(wx.EVT_MENU, self.on_convert_to_nm, m_convert_to_nm)
        self.spectrum_function_items.append(m_convert_to_nm)
        m_resample_spectrum = menu_edit.Append(-1, "Resample spectrum", "Resample wavelength grid")
        self.Bind(wx.EVT_MENU, self.on_resample_spectrum, m_resample_spectrum)
        self.spectrum_function_items.append(m_resample_spectrum)
        m_combine_spectra = menu_edit.Append(-1, "Combine all spectra", "Combine all spectra into one")
        self.Bind(wx.EVT_MENU, self.on_combine_spectra, m_combine_spectra)
        self.spectrum_function_items.append(m_combine_spectra)

        if "generate_spectrum" in dir(sve):
            menu_edit.AppendSeparator()
            m_synthesize = menu_edit.Append(-1, "&Synthesize spectrum", "Synthesize spectrum")
            self.Bind(wx.EVT_MENU, self.on_synthesize, m_synthesize)

        menu_help = wx.Menu()
        m_about = menu_help.Append(-1, "&About\tF1", "About the visual editor")
        self.Bind(wx.EVT_MENU, self.on_about, m_about)

        self.menubar.Append(menu_file, "&File")
        self.menubar.Append(menu_edit, "&Edit")
        self.menubar.Append(menu_help, "&Help")
        self.SetMenuBar(self.menubar)
        self.update_menu_active_spectrum()

    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = wx.Panel(self)

        # Create the mpl Figure and FigCanvas objects.
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 5.0), dpi=self.dpi)
        self.canvas = FigCanvas(self.panel, -1, self.fig)

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

        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

        text_action = wx.StaticText(self.panel, -1, "Action: ", style=wx.ALIGN_LEFT)
        self.radio_button_stats = wx.RadioButton(self.panel, -1, 'Stats', style=wx.RB_GROUP)
        self.radio_button_create = wx.RadioButton(self.panel, -1, 'Create')
        self.radio_button_modify = wx.RadioButton(self.panel, -1, 'Modify')
        self.radio_button_remove = wx.RadioButton(self.panel, -1, 'Remove')
        self.radio_button_stats.SetValue(True)
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_create.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_modify.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_remove.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_action_change, id=self.radio_button_stats.GetId())

        text_element = wx.StaticText(self.panel, -1, "Elements: ", style=wx.ALIGN_LEFT)
        self.radio_button_continuum = wx.RadioButton(self.panel, -1, 'Continuum', style=wx.RB_GROUP)
        self.radio_button_lines = wx.RadioButton(self.panel, -1, 'Lines')
        self.radio_button_linemarks = wx.RadioButton(self.panel, -1, 'Line marks')
        self.radio_button_segments = wx.RadioButton(self.panel, -1, 'Segments')
        self.radio_button_continuum.SetValue(True)
        self.radio_button_linemarks.Enable(False) # Because "Stats" is selected as an action by default
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_continuum.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_lines.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_linemarks.GetId())
        self.Bind(wx.EVT_RADIOBUTTON, self.on_element_change, id=self.radio_button_segments.GetId())

        # Progress bar
        self.gauge = wx.Gauge(self.panel, range=100, size=(100, 20), style=wx.GA_HORIZONTAL)

        self.stats = wx.ListCtrl(self.panel, -1, style=wx.LC_REPORT)
        self.stats.InsertColumn(0, 'Property')
        self.stats.InsertColumn(1, 'Value')
        self.stats.SetColumnWidth(0, 300)
        self.stats.SetColumnWidth(1, 300)

        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas)

        #
        # Layout with box sizers
        #
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox.Add(self.toolbar, 0, wx.EXPAND)
        self.vbox.AddSpacer(10)

        self.hbox = wx.BoxSizer(wx.HORIZONTAL)
        flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL

        self.hbox.Add(text_action, 0, border=3, flag=flags)

        self.hbox.Add(self.radio_button_stats, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_create, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_modify, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_remove, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)

        self.hbox.Add(text_element, 0, border=3, flag=flags)

        self.hbox.Add(self.radio_button_continuum, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_lines, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_linemarks, 0, border=3, flag=flags)
        self.hbox.Add(self.radio_button_segments, 0, border=3, flag=flags)
        self.hbox.AddSpacer(30)

        self.parent_hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.parent_hbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP)

        # Progress bar
        self.parent_hbox.Add(self.gauge, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)

        self.vbox.Add(self.parent_hbox, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP)

        self.vbox.Add(self.stats, 0, flag = wx.ALIGN_LEFT | wx.ALL | wx.TOP | wx.EXPAND)

        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)



    def create_status_bar(self):
        self.statusbar = self.CreateStatusBar()

    def draw_active_spectrum(self):
        if self.active_spectrum != None:
            # Remove spectrum plot if exists
            if self.active_spectrum.plot_id != None:
                self.axes.lines.remove(self.active_spectrum.plot_id)

            # zorder = 1, always in the background
            self.active_spectrum.plot_id = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'], lw=1, color=self.active_spectrum.color, linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor=self.active_spectrum.color, zorder=1, label="[A] "+self.active_spectrum.name)[0]

            # Draw errors
            if self.m_show_errors.IsChecked():
                if self.active_spectrum.errors_plot_id1 != None:
                    self.axes.lines.remove(self.active_spectrum.errors_plot_id1)
                if self.active_spectrum.errors_plot_id2 != None:
                    self.axes.lines.remove(self.active_spectrum.errors_plot_id2)

                # zorder = 1, always in the background
                self.active_spectrum.errors_plot_id1 = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'] + self.active_spectrum.data['err'], lw=1, color=self.active_spectrum.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]
                self.active_spectrum.errors_plot_id2 = self.axes.plot(self.active_spectrum.data['waveobs'], self.active_spectrum.data['flux'] - self.active_spectrum.data['err'], lw=1, color=self.active_spectrum.color, linestyle='-.', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]

        # Put a legend to the right of the current axis
        self.update_legend()

    def remove_drawn_continuum_spectrum(self):
        if self.active_spectrum != None and self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
            self.active_spectrum.continuum_plot_id = None
            self.canvas.draw()

    def remove_continuum_spectrum(self):
        self.remove_drawn_continuum_spectrum()
        self.active_spectrum.continuum_model = None
        self.active_spectrum.continuum_data = None

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

    def draw_continuum_spectrum(self):
        if self.active_spectrum == None:
            return

        # Remove continuum plot if exists
        if self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)

        # zorder = 1, always in the background
        if self.active_spectrum.continuum_data != None:
            self.active_spectrum.continuum_plot_id = self.axes.plot(self.active_spectrum.continuum_data['waveobs'], self.active_spectrum.continuum_data['flux'], lw=1, color='green', linestyle='-', marker='', markersize=1, markeredgewidth=0, markerfacecolor='b', zorder=1)[0]


    def show_marks(self):
        elements = "lines"
        for r in self.region_widgets[elements]:
            r.show_mark()

    def hide_marks(self):
        elements = "lines"
        for r in self.region_widgets[elements]:
            r.hide_mark()

    # Draw all elements in regions array and creates widgets in regions widget array
    # Elements cann be "continuum", "lines" or "segments"
    def draw_regions(self, elements):
        for r in self.region_widgets[elements]:
            r.disconnect_and_remove()
        self.region_widgets[elements] = []

        if elements == "continuum":
            color = 'green'
        elif elements == "lines":
            color = 'yellow'
        else:
            color = 'grey'

        for r in self.regions[elements]:
            axvs = self.axes.axvspan(r['wave_base'], r['wave_top'], facecolor=color, alpha=0.30)
            axvs.zorder = 3

            # Segments always in the background but above spectra
            if elements == "segments":
                axvs.zorder = 2

            if elements == "lines":
                # http://matplotlib.sourceforge.net/examples/pylab_examples/annotation_demo.html
                # http://matplotlib.sourceforge.net/examples/pylab_examples/annotation_demo2.html
                if r['note'] != None and r['note'] != "":
                    note = self.axes.annotate(r['note'], xy=(r['wave_peak'], 1),  xycoords=("data", 'axes fraction'),
                        xytext=(-10, 20), textcoords='offset points',
                        size=8,
                        bbox=dict(boxstyle="round", fc="0.8"),
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                        edgecolor='black'),
                        horizontalalignment='right', verticalalignment='top',
                        annotation_clip=True,
                        )
                else:
                    note = None
                axvline = self.axes.axvline(x = r['wave_peak'], linewidth=1, color='orange')
                region = CustomizableRegion(self, "lines", axvs, mark=axvline, note=note)
            else:
                region = CustomizableRegion(self, elements, axvs)
            region.connect()
            self.region_widgets[elements].append(region)

    def draw_figure_for_first_time(self):
        """ Redraws the figure
        """
        self.axes.grid(True, which="both")
        self.axes.set_title("Spectra", fontsize="10")
        self.axes.set_xlabel("wavelength", fontsize="10")
        self.axes.set_ylabel("flux", fontsize="10")

        for spec in self.spectra:
            # Remove "[A]  " from spectrum name (legend) if it exists
            if self.active_spectrum != None and self.active_spectrum.plot_id != None:
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
            self.active_spectrum = spec
            self.draw_active_spectrum()
        self.draw_regions("continuum")
        self.draw_regions("lines")
        self.draw_regions("segments")
        self.canvas.draw()

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

    def update_title(self):
        title = self.title

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
        self.SetTitle(title)

    # Change saved status and title for the concrete element
    def regions_changed(self, element_type):
        self.not_saved[element_type] = True
        self.update_title()

    # Change saved status and title for the concrete element
    def regions_saved(self, element_type):
        self.not_saved[element_type] = False
        self.update_title()


    def on_close(self, event):
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
                self.tbicon.Destroy()
                self.Destroy()
        else:
            self.tbicon.Destroy()
            self.Destroy()

        # Embeded ipython terminal
        if self.ipython_thread != None:
            self.ipython_thread.join()


    def on_action_change(self, event):
        self.enable_elements()
        self.stats.DeleteAllItems()

        if self.radio_button_create.GetId() == event.GetId():
            self.action = "Create"
        elif self.radio_button_modify.GetId() == event.GetId():
            self.action = "Modify"
        elif self.radio_button_remove.GetId() == event.GetId():
            self.action = "Remove"
        else:
            self.action = "Stats"
            self.radio_button_linemarks.Enable(False)
            if self.radio_button_linemarks.GetValue():
                self.radio_button_continuum.SetValue(True)
                self.elements = "continuum"

    def disable_elements(self):
        self.radio_button_continuum.Enable(False)
        self.radio_button_lines.Enable(False)
        self.radio_button_linemarks.Enable(False)
        self.radio_button_segments.Enable(False)

    def enable_elements(self):
        self.radio_button_continuum.Enable(True)
        self.radio_button_lines.Enable(True)
        self.radio_button_linemarks.Enable(True)
        self.radio_button_segments.Enable(True)

    def on_element_change(self, event):
        if self.radio_button_lines.GetId() == event.GetId():
            self.elements = "lines"
            self.subelements = None
        elif self.radio_button_linemarks.GetId() == event.GetId():
            self.elements = "lines"
            self.subelements = "marks"
        elif self.radio_button_segments.GetId() == event.GetId():
            self.elements = "segments"
            self.subelements = None
        else:
            self.elements = "continuum"
            self.subelements = None

    def on_motion(self, event):
        # Validate it is not in PAN or ZOOM mode
        if event.inaxes == None: return
        if event.inaxes.get_navigate_mode() != None: return
        if self.timeroff.IsRunning(): return
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
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='green', alpha=0.30)
                axvs.zorder = 3
                region = CustomizableRegion(self, "continuum", axvs)
                region.connect()
                self.region_widgets['continuum'].append(region)
            elif self.elements == "lines":
                axvs = self.axes.axvspan(event.xdata - new_halfwidth, event.xdata + new_halfwidth, facecolor='yellow', alpha=0.30)
                axvs.zorder = 3
                axvline = self.axes.axvline(x = event.xdata, linewidth=1, color='orange')

                dlg = wx.TextEntryDialog(self, 'Note for the new line region:','Note')
                dlg.SetValue("")
                if dlg.ShowModal() == wx.ID_OK:
                    note_text = dlg.GetValue()
                    if note_text != "":
                        note = self.axes.annotate(note_text, xy=(event.xdata, 1),  xycoords=("data", 'axes fraction'),
                            xytext=(-10, 20), textcoords='offset points',
                            size=8,
                            bbox=dict(boxstyle="round", fc="0.8"),
                            arrowprops=dict(arrowstyle="->",
                                            connectionstyle="angle,angleA=0,angleB=90,rad=10",
                                            edgecolor='black'),
                            horizontalalignment='right', verticalalignment='top',
                            annotation_clip=True,
                            )
                    else:
                        note = None
                else:
                    note = None
                dlg.Destroy()

                region = CustomizableRegion(self, "lines", axvs, mark=axvline, note=note)
                region.connect()
                self.region_widgets['lines'].append(region)
                self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - new_halfwidth) + " to " + "%.4f" % (event.xdata + new_halfwidth))
            elif self.elements == "segments":
                factor = 100
                axvs = self.axes.axvspan(event.xdata - factor*new_halfwidth, event.xdata + factor*new_halfwidth, facecolor='grey', alpha=0.30)
                # Segments always in the background but above spectra
                axvs.zorder = 2
                region = CustomizableRegion(self, "segments", axvs)
                region.connect()
                self.region_widgets['segments'].append(region)
                self.flash_status_message("Create new region from " + "%.4f" % (event.xdata - factor*new_halfwidth) + " to " + "%.4f" % (event.xdata + factor*new_halfwidth))
            self.canvas.draw()
            self.regions_changed(self.elements)
            ## Automatically change from create to modify, because this will be
            ## the general will of the user
            #~ self.radio_button_create.SetValue(False)
            #~ self.radio_button_modify.SetValue(True)
            #~ self.action = "Modify"

    def update_progress(self, value):
        wx.CallAfter(self.gauge.SetValue, pos=value)

    def update_menu_active_spectrum(self):
        # Remove everything
        for i in np.arange(self.menu_active_spectrum.MenuItemCount):
            self.menu_active_spectrum.Remove(i+1)

        # Add as many options as spectra
        for i in np.arange(len(self.spectra)):
            if self.spectra[i] == None:
                continue

            spec_element = self.menu_active_spectrum.Append(i+1, self.spectra[i].name, kind=wx.ITEM_RADIO)
            self.Bind(wx.EVT_MENU, self.on_change_active_spectrum, spec_element)

            if self.active_spectrum == self.spectra[i]:
                spec_element.Check()

        if len(self.spectra) > 0:
            for item in self.spectrum_function_items:
                item.Enable(True)
        else:
            for item in self.spectrum_function_items:
                item.Enable(False)

    def update_legend(self):
        if len(self.spectra) == 0:
            self.axes.legend_ = None
        else:
            leg = self.axes.legend(loc='upper right', bbox_to_anchor=(1.4, 1.1), ncol=1, shadow=False)
            ltext  = leg.get_texts()
            plt.setp(ltext, fontsize='8')

    def add_stats(self, k, v):
        num_items = self.stats.GetItemCount()
        self.stats.InsertStringItem(num_items, k)
        self.stats.SetStringItem(num_items, 1, str(v))

    def update_stats(self, region):
        if not self.check_active_spectrum_exists():
            return

        self.stats.DeleteAllItems()

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
            rms = np.mean(residuals) + np.std(residuals)
            self.add_stats("Continuum fit root mean square (RMS)", "%.4f" % rms)



    def open_file(self, elements):
        file_choices = "All|*"

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
            if elements == "spectra":
                style = wx.OPEN|wx.FD_MULTIPLE
            else:
                style = wx.OPEN
            dlg = wx.FileDialog(
                self,
                message="Open %s..." % elements,
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=style)

            if dlg.ShowModal() == wx.ID_OK:

                try:
                    if elements == "spectra":
                        paths = dlg.GetPaths()
                        filenames = dlg.GetFilenames()
                        some_does_not_exists = False
                        for i, path in enumerate(paths):
                            if not os.path.exists(path):
                                msg = 'File %s does not exist.' % filenames[i]
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
                            self.update_menu_active_spectrum()
                            self.draw_active_spectrum()
                        self.update_scale()

                        if len(paths) == 1:
                            self.flash_status_message("Opened file %s" % paths[0])
                        else:
                            self.flash_status_message("Opened %i spectra files" % len(paths))
                    else:
                        path = dlg.GetPath()
                        if not os.path.exists(path):
                            msg = 'File %s does not exist.' % dlg.GetFilename()
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

    def on_open_spectra(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("spectra")

    def on_open_continuum(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("continuum")

    def on_open_lines(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("lines")

    def on_open_segments(self, event):
        if self.check_operation_in_progress():
            return
        self.open_file("segments")

    def on_close_spectrum(self, event):
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
        if self.active_spectrum != None and self.active_spectrum.errors_plot_id1 != None:
            self.axes.lines.remove(self.active_spectrum.errors_plot_id1)
            self.active_spectrum.errors_plot_id1 = None
        if self.active_spectrum != None and self.active_spectrum.errors_plot_id2 != None:
            self.axes.lines.remove(self.active_spectrum.errors_plot_id2)
            self.active_spectrum.errors_plot_id2 = None
        # Remove fitted continuum if it exists
        if self.active_spectrum != None and self.active_spectrum.continuum_plot_id != None:
            self.axes.lines.remove(self.active_spectrum.continuum_plot_id)
            self.active_spectrum.continuum_plot_id = None
            self.active_spectrum.continuum_model = None
            self.active_spectrum.continuum_data = None
        # Remove fitted lines if they exist
        for region in self.region_widgets["lines"]:
            if region.line_model.has_key(self.active_spectrum):
                if region.line_plot_id[self.active_spectrum] != None:
                    self.axes.lines.remove(region.line_plot_id[self.active_spectrum])
                del region.line_plot_id[self.active_spectrum]
                del region.line_model[self.active_spectrum]
                del region.line_extra[self.active_spectrum]
        if len(self.spectra) == 0:
            self.active_spectrum = None
        else:
            self.active_spectrum = self.spectra[0]

        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.update_scale()
        self.flash_status_message("Spectrum closed.")

    def on_close_all_spectra(self, event):
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

        self.update_menu_active_spectrum()
        self.update_title()
        self.draw_active_spectrum()
        self.draw_continuum_spectrum()
        self.draw_fitted_lines()
        self.update_scale()
        self.flash_status_message("All spectra closed.")

    def on_save_plot(self, event):
        if self.check_operation_in_progress():
            return
        file_choices = "PNG (*.png)|*.png"

        elements = "spectra"
        if self.filenames[elements] != None:
            filename = self.filenames[elements].split('/')[-1]
            filename_length = len(filename)
            dirname = self.filenames[elements][:-filename_length]
            filename = filename + ".png"
        else:
            filename = elements + ".png"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self,
                message="Save plot as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)

            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    msg = 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename()
                    title = 'File already exists'
                    if not self.question(title, msg):
                        continue # Give the oportunity to select a new file name
                self.canvas.print_figure(path, dpi=self.dpi)
                self.flash_status_message("Saved to %s" % path)
                action_ended = True
            else:
                self.flash_status_message("Not saved.")
                action_ended = True


    def on_save_spectrum(self, event):
        if self.check_operation_in_progress():
            return
        file_choices = "All|*"

        if not self.check_active_spectrum_exists():
            return

        if self.active_spectrum.path != None:
            filename = self.active_spectrum.path.split('/')[-1]
            filename_length = len(filename)
            dirname = self.active_spectrum.path[:-filename_length]
        else:
            filename = self.active_spectrum.name + ".txt"
            dirname = os.getcwd()

        action_ended = False
        while not action_ended:
            dlg = wx.FileDialog(
                self,
                message="Save spectrum as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)

            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    msg = 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename()
                    title = 'File already exists'
                    if not self.question(title, msg):
                        continue # Give the oportunity to select a new file name

                self.status_message("Saving %s..." % path)
                # Save, compress if the filename ends with ".gz"
                sve.write_spectrum(self.active_spectrum.data, path, compress=(path.split('.')[-1] == "gz"))
                self.active_spectrum.not_saved = False
                self.update_title()

                # Change name and path
                self.active_spectrum.path = path
                self.active_spectrum.name = path.split('/')[-1]
                self.active_spectrum.plot_id.set_label("[A] " + self.active_spectrum.name)
                self.update_legend()
                self.canvas.draw()

                # Menu active spectrum
                for item in self.menu_active_spectrum.GetMenuItems():
                    if item.IsChecked():
                        item.SetItemLabel(self.active_spectrum.name)
                        break

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
            dlg = wx.FileDialog(
                self,
                message="Save regions as...",
                defaultDir=dirname,
                defaultFile=filename,
                wildcard=file_choices,
                style=wx.SAVE)

            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                if os.path.exists(path):
                    msg = 'Are you sure you want to overwrite the file %s?' % dlg.GetFilename()
                    title = 'File already exists'
                    if not self.question(title, msg):
                        continue # Give the oportunity to select a new file name
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


    def on_save_continuum_regions(self, event):
        if self.check_operation_in_progress():
            return
        self.save_regions("continuum")

    def on_save_line_regions(self, event):
        if self.check_operation_in_progress():
            return
        self.save_regions("lines")

    def on_save_segments(self, event):
        if self.check_operation_in_progress():
            return
        self.save_regions("segments")

    def text2float(self, value, msg):
        result = None
        try:
            result = float(value)
        except ValueError as e:
            title = 'Bad value format'
            self.error(title, msg)

        return result

    def on_remove_fitted_lines(self, event):
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

    def on_remove_continuum_regions(self, event):
        if self.check_operation_in_progress():
            return
        elements = "continuum"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned continuum regions")


    def on_remove_line_masks(self, event):
        if self.check_operation_in_progress():
            return
        elements = "lines"
        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements)
        self.flash_status_message("Cleaned line masks")

    def on_remove_segments(self, event):
        if self.check_operation_in_progress():
            return
        elements = "segments"
        self.remove_regions(elements)
        self.flash_status_message("Cleaned segments")

    def on_degrade_resolution(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        # Give priority to the resolution determined by the velocity profile relative to telluric lines
        if self.active_spectrum.resolution_telluric == 0.0:
            R = self.active_spectrum.resolution_atomic
        else:
            R = self.active_spectrum.resolution_telluric
        dlg = DegradeResolutionDialog(self, -1, "Degrade spectrum resolution", R, R/2)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        from_resolution = self.text2float(dlg.from_resolution.GetValue(), 'Initial resolution value is not a valid one.')
        to_resolution = self.text2float(dlg.to_resolution.GetValue(), 'Final resolution value is not a valid one.')
        dlg.Destroy()

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
        wx.CallAfter(self.status_message, "Degrading spectrum resolution...")
        if from_resolution == 0:
            # Smooth
            convolved_spectrum = sve.convolve_spectrum(self.active_spectrum.data, to_resolution, frame=self)
        else:
            convolved_spectrum = sve.convolve_spectrum(self.active_spectrum.data, from_resolution, to_resolution=to_resolution, frame=self)
        wx.CallAfter(self.on_degrade_resolution_finnish, convolved_spectrum, from_resolution, to_resolution)

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

    def on_clean_spectrum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        dlg = CleanSpectrumDialog(self, -1, "Clean fluxes and errors", 0.0, 1.2, 0.0, 1.2)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        flux_base = self.text2float(dlg.flux_base.GetValue(), 'Base flux value is not a valid one.')
        flux_top = self.text2float(dlg.flux_top.GetValue(), 'Top flux value is not a valid one.')
        err_base = self.text2float(dlg.err_base.GetValue(), 'Base error value is not a valid one.')
        err_top = self.text2float(dlg.err_top.GetValue(), 'Top error value is not a valid one.')
        dlg.Destroy()

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

        ffilter = (self.active_spectrum.data['flux'] > flux_base) & (self.active_spectrum.data['flux'] <= flux_top)
        efilter = (self.active_spectrum.data['err'] > err_base) & (self.active_spectrum.data['err'] <= err_top)
        wfilter = np.logical_and(ffilter, efilter)
        self.active_spectrum.data = self.active_spectrum.data[wfilter]
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum cleaned.")

    def on_cut_spectrum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        dlg = CutSpectrumDialog(self, -1, "Wavelength range reduction", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2))
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        wave_base = self.text2float(dlg.wave_base.GetValue(), 'Base wavelength value is not a valid one.')
        wave_top = self.text2float(dlg.wave_top.GetValue(), 'Top wavelength value is not a valid one.')
        dlg.Destroy()

        if wave_base == None or wave_top == None or wave_top <= wave_base:
            self.flash_status_message("Bad value.")
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to cut it now anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return

        wfilter = (self.active_spectrum.data['waveobs'] >= wave_base) & (self.active_spectrum.data['waveobs'] <= wave_top)
        if len(self.active_spectrum.data[wfilter]) == 0:
            msg = "This action cannot be done since it would produce a spectrum without measurements."
            title = "Wrong wavelength range"
            self.error(title, msg)
            self.flash_status_message("Bad value.")
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

    def on_resample_spectrum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        dlg = ResampleSpectrumDialog(self, -1, "Resample spectrum", np.round(np.min(self.active_spectrum.data['waveobs']), 2), np.round(np.max(self.active_spectrum.data['waveobs']), 2), 0.001)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        wave_base = self.text2float(dlg.wave_base.GetValue(), 'Base wavelength value is not a valid one.')
        wave_top = self.text2float(dlg.wave_top.GetValue(), 'Top wavelength value is not a valid one.')
        wave_step = self.text2float(dlg.wave_step.GetValue(), 'Wavelength step value is not a valid one.')
        dlg.Destroy()

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

        thread = threading.Thread(target=self.on_resample_spectrum_thread, args=(wave_base, wave_top, wave_step,))
        thread.setDaemon(True)
        thread.start()

    def on_resample_spectrum_thread(self, wave_base, wave_top, wave_step):
        # Homogenize
        wx.CallAfter(self.status_message, "Resampling spectrum...")
        wx.CallAfter(self.update_progress, 10)
        xaxis = np.arange(wave_base, wave_top, wave_step)
        resampled_spectrum_data = sve.resample_spectrum(self.active_spectrum.data, xaxis, linear=True, frame=self)
        self.active_spectrum.data = resampled_spectrum_data
        wx.CallAfter(self.on_resample_spectrum_finnish)

    def on_resample_spectrum_finnish(self):
        self.active_spectrum.not_saved = True

        self.draw_active_spectrum()
        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectrum resampled.")
        self.operation_in_progress = False

    def on_combine_spectra(self, event):
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

        dlg = CombineSpectraDialog(self, -1, "Resample & combine spectrum", np.round(wave_base, 2), np.round(wave_top, 2), 0.001)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        wave_base = self.text2float(dlg.wave_base.GetValue(), 'Base wavelength value is not a valid one.')
        wave_top = self.text2float(dlg.wave_top.GetValue(), 'Top wavelength value is not a valid one.')
        wave_step = self.text2float(dlg.wave_step.GetValue(), 'Wavelength step value is not a valid one.')
        operation_median = dlg.radio_button_median.GetValue()
        operation_mean = dlg.radio_button_mean.GetValue()
        operation_subtract = dlg.radio_button_subtract.GetValue()
        operation_add = dlg.radio_button_add.GetValue()
        operation_divide = dlg.radio_button_divide.GetValue()
        dlg.Destroy()

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
            wx.CallAfter(self.status_message, "Resampling spectrum %i of %i (%s)..." % (i+1, total, spec.name))
            if spec == self.active_spectrum:
                active = i
            resampled_spectrum_data = sve.resample_spectrum(spec.data, xaxis, frame=self)
            resampled_spectra.append(resampled_spectrum_data)
            i += 1

        # Combine
        wx.CallAfter(self.status_message, "Combining spectra...")
        wx.CallAfter(self.update_progress, 10)

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

        wx.CallAfter(self.on_combine_spectra_finnish, combined_spectrum, combined_spectrum_name)

    def on_combine_spectra_finnish(self, combined_spectrum, combined_spectrum_name):
        # Remove "[A]  " from spectrum name (legend)
        self.active_spectrum.plot_id.set_label(self.active_spectrum.name)

        # Add plot
        name = self.get_name(combined_spectrum_name) # If it already exists, add a suffix
        color = self.get_color()
        self.active_spectrum = Spectrum(combined_spectrum, name, color=color)
        self.active_spectrum.not_saved = True
        self.spectra.append(self.active_spectrum)
        self.draw_active_spectrum()

        self.update_menu_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Spectra combined.")
        self.operation_in_progress = False


    def on_continuum_normalization(self, event):
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

    def on_fit_lines(self, event):
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

        vel_telluric = self.active_spectrum.velocity_telluric
        vel_atomic = self.active_spectrum.velocity_atomic
        dlg = FitLinesDialog(self, -1, "Fit lines and cross-match with VALD linelist", vel_atomic, vel_telluric)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        vel_telluric = self.text2float(dlg.vel_telluric.GetValue(), 'Velocity respect to telluric lines is not a valid one.')
        vel_atomic = self.text2float(dlg.vel_atomic.GetValue(), 'Velocity respect to atomic lines is not a valid one.')
        dlg.Destroy()
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
        self.update_numpy_arrays_from_widgets("lines")
        vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")
        telluric_linelist_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
        linemasks = sve.fit_lines(self.regions["lines"], self.active_spectrum.data, self.active_spectrum.continuum_model, vel_atomic, vel_telluric, vald_linelist_file, chemical_elements_file, molecules_file, telluric_linelist_file, discard_gaussian=False, discard_voigt=True, frame=self)
        wx.CallAfter(self.on_fit_lines_finnish, linemasks)

    def on_fit_lines_finnish(self, linemasks):
        elements = "lines"
        self.remove_fitted_lines() # If they exist
        self.remove_regions(elements, check_not_saved=False)

        if linemasks != None and len(linemasks) > 0:
            total_regions = len(linemasks)
            line_regions = np.recarray((total_regions, ), dtype=[('wave_peak', float),('wave_base', float), ('wave_top', float), ('note', '|S100')])
            line_regions['wave_peak'] = linemasks['mu']
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

            self.regions[elements] = line_regions
            self.draw_regions(elements)

            # Fitted Gaussian lines
            i = 0
            for region in self.region_widgets["lines"]:
                region.line_model[self.active_spectrum] = line_models[i]
                region.line_extra[self.active_spectrum] = line_extras[i]
                i += 1
            self.draw_fitted_lines()

            self.not_saved[elements] = True
            self.update_title()
            self.flash_status_message("%i line masks found/fitted!" % (len(line_regions)))
        else:
            self.flash_status_message("No line masks found/fitted!")
        self.canvas.draw()
        self.operation_in_progress = False


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

    def on_change_active_spectrum(self, event):
        if self.check_operation_in_progress():
            return

        # Remove current continuum from plot if exists
        self.remove_drawn_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        self.remove_drawn_fitted_lines()

        # Determine the new spectrum and draw its continuum if exists
        i = 0
        for item in self.menu_active_spectrum.GetMenuItems():
            if item.IsChecked():
                # Remove "[A]  " from spectrum name (legend)
                self.active_spectrum.plot_id.set_label(self.active_spectrum.name)
                # Change active spectrum
                self.active_spectrum = self.spectra[i]
                self.draw_active_spectrum()

                self.draw_continuum_spectrum()
                self.draw_fitted_lines()
                self.canvas.draw()
                self.flash_status_message("Active spectrum: %s" % self.active_spectrum.name)
                break
            i += 1

    def on_remove_fitted_continuum(self, event):
        if self.check_operation_in_progress():
            return
        self.remove_continuum_spectrum()

    def on_fit_continuum(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        # Initial recommendation: 1 knot every 10 nm
        nknots = np.max([1, int((np.max(self.active_spectrum.data['waveobs']) - np.min(self.active_spectrum.data['waveobs'])) / 10)])
        median_wave_range=0.1
        max_wave_range=1

        dlg = FitContinuumDialog(self, -1, "Properties for fitting continuum", nknots=nknots, median_wave_range=median_wave_range,
max_wave_range=max_wave_range)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return


        nknots = self.text2float(dlg.nknots.GetValue(), 'Number of knots value is not a valid one.')
        median_wave_range = self.text2float(dlg.median_wave_range.GetValue(), 'Median wave range value is not a valid one.')
        max_wave_range = self.text2float(dlg.max_wave_range.GetValue(), 'Max wave range value is not a valid one.')
        in_continuum = dlg.radio_button_continuum.GetValue()
        dlg.Destroy()

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

    def get_spectrum_from_model(self, model, spectrum_wave_grid):
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
                self.update_numpy_arrays_from_widgets("continuum")
                self.active_spectrum.continuum_model = sve.fit_continuum(self.active_spectrum.data, segments=self.regions["continuum"] , nknots=nknots, median_wave_range=median_wave_range, max_wave_range=max_wave_range)
            else:
                self.active_spectrum.continuum_model = sve.fit_continuum(self.active_spectrum.data, nknots=nknots, median_wave_range=median_wave_range, max_wave_range=max_wave_range)
            self.active_spectrum.continuum_data = self.get_spectrum_from_model(self.active_spectrum.continuum_model, self.active_spectrum.data['waveobs'])
            wx.CallAfter(self.on_fit_continuum_finish, nknots)
        except Exception as e:
            self.operation_in_progress = False
            wx.CallAfter(self.flash_status_message, "Not enough data points to fit, reduce the number of nknots or increase the spectrum regions.")


    def on_fit_continuum_finish(self, nknots):
        self.draw_continuum_spectrum()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Continuum fitted with %s knots uniform Spline model." % str(nknots))

    def on_find_continuum(self, event):
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

        dlg = FindContinuumDialog(self, -1, "Properties for finding continuum regions", self.find_continuum_regions_wave_step, self.find_continuum_regions_sigma, self.find_continuum_regions_max_continuum_diff)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        #~ resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
        resolution = None
        fixed_wave_step = self.text2float(dlg.fixed_wave_step.GetValue(), 'Maximum standard deviation value is not a valid one.')
        sigma = self.text2float(dlg.sigma.GetValue(), 'Maximum standard deviation value is not a valid one.')
        max_continuum_diff = self.text2float(dlg.max_continuum_diff.GetValue(), 'Maximum continuum difference value is not a valid one.')
        in_segments = dlg.radio_button_segments.GetValue()
        dlg.Destroy()

        if fixed_wave_step == None or sigma == None or max_continuum_diff == None:
            wx.CallAfter(self.flash_status_message, "Bad value.")
            return
        # Save values
        self.find_continuum_regions_wave_step = fixed_wave_step
        self.find_continuum_regions_sigma = sigma
        self.find_continuum_regions_max_continuum_diff = max_continuum_diff
        # Convert from % to over 1
        max_continuum_diff = max_continuum_diff / 100

        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            wx.CallAfter(self.flash_status_message, "No segments found.")
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
        wx.CallAfter(self.status_message, "Finding continuum regions...")
        if in_segments:
            self.update_numpy_arrays_from_widgets("segments")
            continuum_regions = sve.find_continuum(self.active_spectrum.data, resolution, segments=self.regions["segments"], max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)
        else:
            continuum_regions = sve.find_continuum(self.active_spectrum.data, resolution, max_std_continuum = sigma, continuum_model = self.active_spectrum.continuum_model, max_continuum_diff=max_continuum_diff, fixed_wave_step=fixed_wave_step, frame=self)

        wx.CallAfter(self.on_find_continuum_finish, continuum_regions)

    def on_find_continuum_finish(self, continuum_regions):
        elements = "continuum"
        self.regions[elements] = continuum_regions
        self.draw_regions(elements)
        self.not_saved[elements] = True
        self.update_title()
        self.canvas.draw()
        self.operation_in_progress = False
        self.flash_status_message("Automatic finding of continuum regions ended.")

    def on_find_lines(self, event):
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
        dlg = FindLinesDialog(self, -1, "Properties for finding line masks", self.find_lines_min_depth, self.find_lines_max_depth, vel_atomic=vel_atomic, vel_telluric=vel_telluric, resolution=R, elements="Fe 1, Fe 2")
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        max_depth = self.text2float(dlg.max_depth.GetValue(), 'Maximum depth value is not a valid one.')
        min_depth = self.text2float(dlg.min_depth.GetValue(), 'Minimum depth value is not a valid one.')
        elements = dlg.elements.GetValue()
        resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
        vel_atomic = self.text2float(dlg.vel_atomic.GetValue(), 'Velocity respect to atomic lines is not a valid one.')
        vel_telluric = self.text2float(dlg.vel_telluric.GetValue(), 'Velocity respect to telluric lines is not a valid one.')
        discard_tellurics = dlg.discard_tellurics.IsChecked()
        in_segments = dlg.radio_button_segments.GetValue()
        dlg.Destroy()

        if max_depth == None or min_depth == None or resolution == None or vel_atomic == None or vel_telluric == None or max_depth <= min_depth or max_depth <= 0 or min_depth < 0 or resolution <= 0:
            wx.CallAfter(self.flash_status_message, "Bad value.")
            return
        ## Save values
        self.find_lines_max_depth = max_depth
        self.find_lines_min_depth = min_depth
        self.active_spectrum.velocity_atomic = vel_atomic
        self.active_spectrum.velocity_telluric = vel_telluric
        self.active_spectrum.resolution_telluric = resolution

        if in_segments and (self.region_widgets["segments"] == None or len(self.region_widgets["segments"]) == 0):
            wx.CallAfter(self.flash_status_message, "No segments found.")
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
        wx.CallAfter(self.status_message, "Smoothing spectrum...")
        smoothed_spectrum = sve.convolve_spectrum(spectrum, 2*resolution, frame=self)


        wx.CallAfter(self.status_message, "Generating line masks, fitting gaussians and matching VALD lines...")
        logging.info("Generating line masks, fitting gaussians and matching VALD lines...")
        vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
        chemical_elements_file = resource_path("input/abundances/chemical_elements_symbols.dat")
        molecules_file = resource_path("input/abundances/molecular_symbols.dat")
        telluric_linelist_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
        linemasks = sve.find_linemasks(spectrum, self.active_spectrum.continuum_model, vald_linelist_file, chemical_elements_file, molecules_file, telluric_linelist_file, minimum_depth=min_depth, maximum_depth=max_depth, smoothed_spectrum=smoothed_spectrum, discard_gaussian = False, discard_voigt = True, vel_atomic=vel_atomic, vel_telluric=vel_telluric, frame=self)

        # If no peaks found, just finnish
        if linemasks == None or len(linemasks) == 0:
            wx.CallAfter(self.on_find_lines_finish, None)
            return

        logging.info("Applying filters to discard bad line masks...")
        wx.CallAfter(self.status_message, "Applying filters to discard bad line masks...")

        #rejected_by_atomic_line_not_found = (linemasks['VALD_wave_peak'] == 0)
        rejected_by_telluric_line = (linemasks['telluric_wave_peak'] != 0)

        discarded = linemasks['wave_peak'] < 0 # All to false

        # In case it is specified, select only given elements
        if elements != "":
            elements = elements.split(",")
            select = linemasks['element'] == "NONE" # All to false
            for element in elements:
                select = np.logical_or(select, linemasks['element'] == element.strip().capitalize())
            discarded = np.logical_or(discarded, np.logical_not(select))

        if discard_tellurics:
            discarded = np.logical_or(discarded, rejected_by_telluric_line)
        #discarded = np.logical_or(discarded, rejected_by_atomic_line_not_found)
        if in_segments:
            # Identify linemasks with too big wavelength range (probably star in between segments)
            wave_diff = (linemasks['wave_top'] - linemasks['wave_base']) / (linemasks['top'] - linemasks['base'])
            wave_diff_selected, wave_diff_selected_filter = sve.sigma_clipping(wave_diff[~discarded], sig=3, meanfunc=np.median) # Discard outliers
            accepted_index = np.arange(len(wave_diff))[~discarded][wave_diff_selected_filter]
            rejected_by_wave_gaps = wave_diff < 0 # Create an array of booleans
            rejected_by_wave_gaps[:] = True # Initialize
            rejected_by_wave_gaps[accepted_index] = False
            discarded = np.logical_or(discarded, rejected_by_wave_gaps)

        linemasks = linemasks[~discarded]
        total_regions = len(linemasks)

        # If no regions found, just finnish
        if total_regions == 0:
            wx.CallAfter(self.on_find_lines_finish, None, None, None)
            return


        wx.CallAfter(self.on_find_lines_finish, linemasks)

    def on_find_lines_finish(self, linemasks):
        self.on_fit_lines_finnish(linemasks)

    def on_determine_barycentric_vel(self, event):
        dlg = DetermineBarycentricCorrectionDialog(self, -1, "Barycentric velocity determination", self.day, self.month, self.year, self.hours, self.minutes, self.seconds, self.ra_hours, self.ra_minutes, self.ra_seconds, self.dec_degrees, self.dec_minutes, self.dec_seconds)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        day = self.text2float(dlg.day.GetValue(), 'Day value is not a valid one.')
        month = self.text2float(dlg.month.GetValue(), 'Month value is not a valid one.')
        year = self.text2float(dlg.year.GetValue(), 'Year value is not a valid one.')

        hours = self.text2float(dlg.hours.GetValue(), 'Hours value is not a valid one.')
        minutes = self.text2float(dlg.minutes.GetValue(), 'Minutes value is not a valid one.')
        seconds = self.text2float(dlg.seconds.GetValue(), 'Seconds value is not a valid one.')

        ra_hours = self.text2float(dlg.ra_hours.GetValue(), 'RA hours value is not a valid one.')
        ra_minutes = self.text2float(dlg.ra_minutes.GetValue(), 'RA minutes value is not a valid one.')
        ra_seconds = self.text2float(dlg.ra_seconds.GetValue(), 'RA seconds value is not a valid one.')

        dec_degrees = self.text2float(dlg.dec_degrees.GetValue(), 'DEC degrees value is not a valid one.')
        dec_minutes = self.text2float(dlg.dec_minutes.GetValue(), 'DEC minutes value is not a valid one.')
        dec_seconds = self.text2float(dlg.dec_seconds.GetValue(), 'DEC seconds value is not a valid one.')

        if None in [day, month, year, hours, minutes, seconds, ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds]:
            self.flash_status_message("Bad value.")
            return

        self.day = day
        self.month = month
        self.year = year
        self.hours = hours
        self.minutes = minutes
        self.seconds = seconds
        self.ra_hours = ra_hours
        self.ra_minutes = ra_minutes
        self.ra_seconds = ra_seconds
        self.dec_degrees = dec_degrees
        self.dec_minutes = dec_minutes
        self.dec_seconds = dec_seconds

        # Project velocity toward star
        self.barycentric_vel = sve.calculate_barycentric_velocity_correction((year, month, day, hours, minutes, seconds), (ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds))

        msg = "Barycentric velocity determined: " + str(self.barycentric_vel) + " km/s"
        title = "Barycentric velocity"
        self.info(title, msg)
        self.flash_status_message(msg)

    def on_estimate_snr(self, event):
        if self.check_operation_in_progress():
            return

        if self.active_spectrum.snr != None:
            msg = "Previous global SNR: %.2f. Re-estimate again? " % self.active_spectrum.snr
            title = "Signal-to-Noise Ratio"
            if not self.question(title, msg):
                return

        dlg = EstimateSNRDialog(self, -1, "Properties for estimating SNR", num_points=10)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        num_points = self.text2float(dlg.num_points.GetValue(), 'Number of points value is not a valid one.')
        estimate_from_flux = dlg.radio_button_from_flux.GetValue()
        dlg.Destroy()

        if num_points == None:
            self.flash_status_message("Bad value.")
            return

        if estimate_from_flux:
            self.operation_in_progress = True
            thread = threading.Thread(target=self.on_estimate_snr_thread, args=(num_points,))
            thread.setDaemon(True)
            thread.start()
        else:
            efilter = self.active_spectrum.data['err'] > 0
            spec = self.active_spectrum.data[efilter]
            if len(spec) > 1:
                estimated_snr = np.mean(spec['flux'] / spec['err'])
                self.on_estimate_snr_finnish(estimated_snr)
            else:
                msg = 'All value errors are set to zero or negative numbers'
                title = "SNR estimation error"
                self.error(title, msg)


    def on_estimate_snr_thread(self, num_points):
        wx.CallAfter(self.status_message, "Estimating SNR for the whole spectrum...")
        estimated_snr = sve.estimate_snr(self.active_spectrum.data['flux'], num_points=num_points, frame=self)
        wx.CallAfter(self.on_estimate_snr_finnish, estimated_snr)

    def on_estimate_snr_finnish(self, estimated_snr):
        self.active_spectrum.snr = estimated_snr
        msg = "Estimated global SNR: %.2f" % self.active_spectrum.snr
        title = "Global Signal-to-Noise Ratio"
        self.info(title, msg)
        self.flash_status_message(msg)
        self.operation_in_progress = False

    def on_estimate_errors(self, event):
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

        if self.active_spectrum.snr != None:
            snr = self.active_spectrum.snr
        else:
            snr = 10.0

        dlg = EstimateErrorsDialog(self, -1, "Estimate spectrum errors", snr)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        snr = self.text2float(dlg.snr.GetValue(), 'SNR value is not a valid one.')
        dlg.Destroy()
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

    def on_show_errors(self, event):
        if self.check_operation_in_progress():
            self.m_show_errors.Check(False)
            return
        if self.m_show_errors.IsChecked():
            self.draw_errors_spectra()
        else:
            self.remove_drawn_errors_spectra()

    def on_determine_velocity_atomic(self, event):
        self.on_determine_velocity(event, relative_to_atomic_data = True)

    def on_determine_velocity_telluric(self, event):
        self.on_determine_velocity(event, relative_to_atomic_data = False)

    def on_determine_velocity(self, event, relative_to_atomic_data = True, show_previous_results=True):
        # Check if velocity has been previously determined and show those results
        if show_previous_results and ((relative_to_atomic_data and self.active_spectrum.velocity_profile_atomic_models != None) or (not relative_to_atomic_data and self.active_spectrum.velocity_profile_telluric_models != None)) :
            if relative_to_atomic_data:
                xcoord = self.active_spectrum.velocity_profile_atomic_xcoord
                fluxes = self.active_spectrum.velocity_profile_atomic_fluxes
                models = self.active_spectrum.velocity_profile_atomic_models
                num_used_lines = self.active_spectrum.velocity_profile_atomic_num_used_lines
                velocity_step = self.active_spectrum.velocity_profile_atomic_rv_step
                telluric_fwhm = 0.0
                snr = self.active_spectrum.velocity_profile_atomic_snr
                title = "Velocity profile relative to atomic lines"
            else:
                xcoord = self.active_spectrum.velocity_profile_telluric_xcoord
                fluxes = self.active_spectrum.velocity_profile_telluric_fluxes
                models = self.active_spectrum.velocity_profile_telluric_models
                num_used_lines = self.active_spectrum.velocity_profile_telluric_num_used_lines
                velocity_step = self.active_spectrum.velocity_profile_telluric_rv_step
                telluric_fwhm = self.active_spectrum.velocity_profile_telluric_fwhm_correction
                snr = 0.0
                title = "Velocity profile relative to telluric lines"

            dlg = VelocityProfileDialog(self, -1, title, xcoord, fluxes, models, num_used_lines, velocity_step, telluric_fwhm, snr=snr)
            dlg.ShowModal()
            recalculate = dlg.recalculate
            dlg.Destroy()
            if not recalculate:
                return

        if self.check_operation_in_progress():
            return

        if relative_to_atomic_data:
            linelist = self.linelist_atomic
            velocity_lower_limit = self.velocity_atomic_lower_limit
            velocity_upper_limit = self.velocity_atomic_upper_limit
            velocity_step = self.velocity_atomic_step
        else:
            linelist = self.linelist_telluric
            velocity_lower_limit = self.velocity_telluric_lower_limit
            velocity_upper_limit = self.velocity_telluric_upper_limit
            velocity_step = self.velocity_telluric_step

        dlg = DetermineVelocityDialog(self, -1, "Velocity determination", rv_lower_limit=velocity_lower_limit, rv_upper_limit=velocity_upper_limit, rv_step=velocity_step)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        rv_lower_limit = self.text2float(dlg.rv_lower_limit.GetValue(), 'Velocity lower limit is not a valid one.')
        rv_upper_limit = self.text2float(dlg.rv_upper_limit.GetValue(), 'Velocity upper limit is not a valid one.')
        rv_step = self.text2float(dlg.rv_step.GetValue(), 'Velocity step is not a valid one.')
        dlg.Destroy()

        if rv_lower_limit == None or rv_upper_limit == None or rv_step == None:
            self.flash_status_message("Bad value.")
            return

        if rv_lower_limit >= rv_upper_limit:
            msg = "Upper velocity limit should be greater than lower limit"
            title = "Velocity error"
            self.error(title, msg)
            return

        if (np.abs(rv_lower_limit) + np.abs(rv_upper_limit)) <= 4*rv_step:
            msg = "Velocity step too small for the established limits"
            title = "Velocity error"
            self.error(title, msg)
            return

        if relative_to_atomic_data:
            self.velocity_atomic_lower_limit = rv_lower_limit
            self.velocity_atomic_upper_limit = rv_upper_limit
            self.velocity_atomic_step = rv_step
        else:
            self.velocity_telluric_lower_limit = rv_lower_limit
            self.velocity_telluric_upper_limit = rv_upper_limit
            self.velocity_telluric_step = rv_step

        self.operation_in_progress = True
        thread = threading.Thread(target=self.on_determine_velocity_thread, args=(relative_to_atomic_data,))
        thread.setDaemon(True)
        thread.start()

    def filter_telluric_lines(self, linelist_telluric, spectrum, velocity_lower_limit, velocity_upper_limit):
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


    def on_determine_velocity_thread(self, relative_to_atomic_data):
        wx.CallAfter(self.status_message, "Determining velocity...")

        if relative_to_atomic_data:
            if self.linelist_atomic == None:
                vald_linelist_file = resource_path("input/linelists/VALD/VALD.300_1100nm_teff_5770.0_logg_4.40.lst")
                self.linelist_atomic = sve.read_VALD_linelist(vald_linelist_file, minimum_depth=0.0)
            linelist = self.linelist_atomic
            velocity_lower_limit = self.velocity_atomic_lower_limit
            velocity_upper_limit = self.velocity_atomic_upper_limit
            velocity_step = self.velocity_atomic_step
        else:
            if self.linelist_telluric == None:
                telluric_lines_file = resource_path("input/linelists/telluric/standard_atm_air_model.lst")
                self.linelist_telluric = sve.read_telluric_linelist(telluric_lines_file, minimum_depth=0.0)
            velocity_lower_limit = self.velocity_telluric_lower_limit
            velocity_upper_limit = self.velocity_telluric_upper_limit
            velocity_step = self.velocity_telluric_step
            linelist = self.filter_telluric_lines(self.linelist_telluric, self.active_spectrum.data, velocity_lower_limit, velocity_upper_limit)

        xcoord, fluxes, num_used_lines = sve.build_velocity_profile(self.active_spectrum.data, linelist, lower_velocity_limit=velocity_lower_limit, upper_velocity_limit=velocity_upper_limit, velocity_step=velocity_step, frame=self)
        wx.CallAfter(self.on_determine_velocity_finish, xcoord, fluxes, relative_to_atomic_data, num_used_lines, linelist)

    def on_determine_velocity_finish(self, xcoord, fluxes, relative_to_atomic_data, num_used_lines, linelist):
        # Modelize
        if relative_to_atomic_data:
            models = sve.modelize_velocity_profile(xcoord, fluxes)
            accept = sve.select_good_velocity_profile_models(models, xcoord, fluxes)
            if len(models[accept]) == 0:
                models = models[:1]
            else:
                models = models[accept]
        else:
            models = sve.modelize_velocity_profile(xcoord, fluxes, only_one_peak=True)

        if len(models) == 0:
            fwhm = 0.0
            telluric_fwhm = 0.0
            R = 0.0
            velocity = 0.0
            snr = 0.0
            self.flash_status_message("Velocity could not be determined!")
        else:
            # Resolving power
            c = 299792458.0 # m/s
            fwhm = models[0].fwhm()[0] # km/s (because xcoord is already velocity)
            if relative_to_atomic_data:
                telluric_fwhm = 0.0
                R = np.int(c/(1000.0*np.round(fwhm, 2)))
                ##### Estimate SNR for the lines and nearby regions
                ##### It is only of interest in case of atomic lines
                # Build the fluxes for the composite models
                first_time = True
                for model in models:
                    if first_time:
                        line_fluxes = model(xcoord)
                        first_time = False
                        continue

                    current_line_fluxes = model(xcoord)
                    wfilter = line_fluxes > current_line_fluxes
                    line_fluxes[wfilter] = current_line_fluxes[wfilter]
                # Substract the line models conserving the base level
                base_level = np.median(fluxes)
                corrected_fluxes = base_level + fluxes - line_fluxes
                # Estimate SNR for the lines and nearby regions
                snr = np.mean(corrected_fluxes) / np.std(corrected_fluxes)
            else:
                # If telluric lines have been used, we can substract its natural FWHM
                # so that we get the real resolution of the instrument (based on the difference in FWHM)
                c = 299792458.0 # m/s
                telluric_fwhm = np.mean((c / (linelist['wave_peak'] / linelist['fwhm'])) / 1000.0) # km/s
                R = np.int(c/(1000.0*np.round(fwhm - telluric_fwhm, 2)))
                snr = 0.0
            # Velocity
            velocity = np.round(models[0].mu(), 2) # km/s
            # A positive velocity indicates the distance between the objects is or was increasing;
            # A negative velocity indicates the distance between the source and observer is or was decreasing.
            self.flash_status_message("Velocity determined: " + str(velocity) + " km/s")
        self.operation_in_progress = False

        if relative_to_atomic_data:
            self.active_spectrum.velocity_atomic = velocity
            # If the spe
            if R != 0:
                self.active_spectrum.resolution_atomic = R
            self.active_spectrum.velocity_profile_atomic_xcoord = xcoord
            self.active_spectrum.velocity_profile_atomic_fluxes = fluxes
            self.active_spectrum.velocity_profile_atomic_models = models
            self.active_spectrum.velocity_profile_atomic_num_used_lines = num_used_lines
            self.active_spectrum.velocity_profile_atomic_rv_step = self.velocity_atomic_step
            self.active_spectrum.velocity_profile_atomic_snr = snr
            velocity_step = self.velocity_atomic_step
            title = "Velocity profile relative to atomic lines"
        else:
            self.active_spectrum.velocity_telluric = velocity
            self.barycentric_vel = velocity # So it will appear in the barycentric vel. correction dialog
            if R != 0:
                self.active_spectrum.resolution_telluric = R
            self.active_spectrum.velocity_profile_telluric_xcoord = xcoord
            self.active_spectrum.velocity_profile_telluric_fluxes = fluxes
            self.active_spectrum.velocity_profile_telluric_models = models
            self.active_spectrum.velocity_profile_telluric_num_used_lines = num_used_lines
            self.active_spectrum.velocity_profile_telluric_rv_step = self.velocity_atomic_step
            self.active_spectrum.velocity_profile_telluric_fwhm_correction = telluric_fwhm
            velocity_step = self.velocity_telluric_step
            title = "Velocity profile relative to telluric lines"

        dlg = VelocityProfileDialog(self, -1, title, xcoord, fluxes, models, num_used_lines, velocity_step, telluric_fwhm=telluric_fwhm, snr=snr)
        dlg.ShowModal()
        recalculate = dlg.recalculate
        dlg.Destroy()
        if recalculate:
            self.on_determine_velocity(None, relative_to_atomic_data=relative_to_atomic_data, show_previous_results=False)


    def on_correct_rv(self, event):
        vel_type = "radial"
        default_vel = self.active_spectrum.velocity_atomic
        self.on_correct_vel(event, vel_type, default_vel)

    def on_correct_barycentric_vel(self, event):
        vel_type = "barycentric"
        default_vel = self.barycentric_vel
        self.on_correct_vel(event, vel_type, default_vel)

    def on_correct_vel(self, event, vel_type, default_vel):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        dlg = CorrectVelocityDialog(self, -1, vel_type, default_vel)
        dlg.ShowModal()

        if not dlg.action_accepted:
            dlg.Destroy()
            return

        barycentric_vel = self.text2float(dlg.rv.GetValue(), vel_type.capitalize() + ' velocity value is not a valid one.')
        in_regions = dlg.radio_button_regions.GetValue()
        dlg.Destroy()

        if barycentric_vel == None:
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
                        self.regions[elements] = sve.correct_velocity_regions(self.regions[elements], barycentric_vel, with_peak=True)
                    else:
                        self.regions[elements] = sve.correct_velocity_regions(self.regions[elements], barycentric_vel)
                    self.draw_regions(elements)
                    self.not_saved[elements] = False
        else:
            if not self.check_active_spectrum_exists():
                return
            self.active_spectrum.data = sve.correct_velocity(self.active_spectrum.data, barycentric_vel)
            self.active_spectrum.not_saved = True
            self.draw_active_spectrum()
        self.update_title()
        self.canvas.draw()
        self.flash_status_message("Applied a " + vel_type + " velocity correction of %s." % barycentric_vel)

    def on_convert_to_nm(self, event):
        if not self.check_active_spectrum_exists():
            return
        if self.check_operation_in_progress():
            return

        # Check if spectrum is saved
        if self.active_spectrum.not_saved:
            msg = "The active spectrum has not been saved, are you sure you want to divide all the wavelength by 10 anyway?"
            title = "Changes not saved"
            if not self.question(title, msg):
                return
        else:
            msg = "This operation is going to divide all the wavelengths of the active spectrum by 10. Are you sure?"
            title = "Confirmation"
            if not self.question(title, msg):
                return

        self.status_message("Converting to nanometers...")

        # Remove current continuum from plot if exists
        self.remove_continuum_spectrum()

        # Remove current drawn fitted lines if they exist
        # IMPORTANT: Before active_spectrum is modified, if not this routine will not work properly
        self.remove_fitted_lines()

        self.active_spectrum.data['waveobs'] = self.active_spectrum.data['waveobs'] / 10
        self.active_spectrum.not_saved = True
        self.draw_active_spectrum()

        self.update_title()
        self.update_scale()
        self.flash_status_message("Wavelength divided by 10.")

    def on_synthesize(self, event):
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
            resolution = 47000
            wave_step = 0.001
            dlg = SyntheticSpectrumDialog(self, -1, "Synthetic spectrum generator", wave_base, wave_top, wave_step, resolution, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff)
            dlg.ShowModal()

            if not dlg.action_accepted:
                dlg.Destroy()
                return

            teff = self.text2float(dlg.teff.GetValue(), 'Effective temperature value is not a valid one.')
            logg = self.text2float(dlg.logg.GetValue(), 'Gravity (log g) value is not a valid one.')
            MH = self.text2float(dlg.MH.GetValue(), 'Metallicity [M/H] value is not a valid one.')
            microturbulence_vel = self.text2float(dlg.microturbulence_vel.GetValue(), 'Microturbulence velocity value is not a valid one.')
            macroturbulence = self.text2float(dlg.macroturbulence.GetValue(), 'Macroturbulence velocity value is not a valid one.')
            vsini = self.text2float(dlg.vsini.GetValue(), 'Rotation value (v sin(i)) is not a valid one.')
            limb_darkening_coeff = self.text2float(dlg.limb_darkening_coeff.GetValue(), 'Limb darkening coefficient is not a valid one.')
            resolution = self.text2float(dlg.resolution.GetValue(), 'Resolution value is not a valid one.')
            wave_base = self.text2float(dlg.wave_base.GetValue(), 'Wavelength min. value is not a valid one.')
            wave_top = self.text2float(dlg.wave_top.GetValue(), 'Wavelength max. value is not a valid one.')
            wave_step = self.text2float(dlg.wave_step.GetValue(), 'Wavelength step value is not a valid one.')

            in_segments = dlg.radio_button_segments.GetValue() # else in spectrum
            dlg.Destroy()

            if teff == None or logg == None or MH == None or microturbulence_vel == None or resolution == None or wave_base == None or wave_top == None:
                self.flash_status_message("Bad value.")
                return


            if self.modeled_layers_pack == None:
                logging.info("Loading modeled atmospheres...")
                self.status_message("Loading modeled atmospheres...")
                self.modeled_layers_pack = sve.load_modeled_layers_pack(resource_path('input/atmospheres/default.modeled_layers_pack.dump'))

            if not sve.valid_atmosphere_target(self.modeled_layers_pack, teff, logg, MH):
                msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] fall out of theatmospheric models."
                title = 'Out of the atmospheric models'
                self.error(title, msg)
                self.flash_status_message("Bad values.")
                return

            # Prepare atmosphere model
            self.status_message("Interpolating atmosphere model...")
            layers = sve.interpolate_atmosphere_layers(self.modeled_layers_pack, teff, logg, MH)
            atm_filename = sve.write_atmosphere(teff, logg, MH, layers)

            # Generate
            if not in_segments:
                if wave_base >= wave_top:
                    msg = "Bad wavelength range definition, maximum value cannot be lower than minimum value."
                    title = 'Wavelength range'
                    self.error(title, msg)
                    self.flash_status_message("Bad values.")
                    return

                waveobs = np.arange(wave_base, wave_top, wave_step)
            else:
                #in_segments
                if len(self.region_widgets["segments"]) == 0:
                    self.flash_status_message("No segments present for synthetic spectrum generation.")
                    return

                # Build wavelength points from regions
                waveobs = None
                for region in self.region_widgets["segments"]:
                    wave_base = region.get_wave_base()
                    wave_top = region.get_wave_top()

                    new_waveobs = np.arange(wave_base, wave_top, wave_step)
                    if waveobs == None:
                        waveobs = new_waveobs
                    else:
                        waveobs = np.hstack((waveobs, new_waveobs))

            # If wavelength out of the linelist file are used, SPECTRUM starts to generate flat spectrum
            if np.min(waveobs) <= 300.0 or np.max(waveobs) >= 1100.0:
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
            thread = threading.Thread(target=self.on_synthesize_thread, args=(waveobs, atm_filename, teff, logg, MH, microturbulence_vel,  macroturbulence, vsini, limb_darkening_coeff, resolution))
            thread.setDaemon(True)
            thread.start()

    def on_synthesize_thread(self, waveobs, atm_filename, teff, logg, MH, microturbulence_vel, macroturbulence, vsini, limb_darkening_coeff, resolution):
        total_points = len(waveobs)

        synth_spectrum = np.recarray((total_points, ), dtype=[('waveobs', float),('flux', float),('err', float)])
        synth_spectrum['waveobs'] = waveobs

        # waveobs is multiplied by 10.0 in order to be converted from nm to armstrongs
        synth_spectrum['flux'] = sve.generate_spectrum(synth_spectrum['waveobs']*10.0, atm_filename, linelist_file = resource_path("input/linelists/SPECTRUM/default.300_1100nm.lst"), abundances_file = resource_path("input/abundances/default.stdatom.dat"), microturbulence_vel = microturbulence_vel, macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, R=resolution, verbose=1, update_progress_func=self.update_progress)


        synth_spectrum.sort(order='waveobs') # Make sure it is ordered by wavelength

        # Remove atmosphere model temporary file
        os.remove(atm_filename)
        wx.CallAfter(self.on_synthesize_finnish, synth_spectrum, teff, logg, MH, microturbulence_vel)

    def on_synthesize_finnish(self, synth_spectrum, teff, logg, MH, microturbulence_vel):
        # Draw new lines
        self.draw_fitted_lines()
        self.operation_in_progress = False
        self.flash_status_message("Lines fitted.")
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
        self.active_spectrum.not_saved = True
        self.update_title()
        self.update_menu_active_spectrum()
        self.draw_active_spectrum()
        self.update_scale()

        self.canvas.draw()

        self.operation_in_progress = False
        self.flash_status_message("Synthetic spectrum generated!")


    def on_exit(self, event):
        self.on_close(event)

    def on_about(self, event):
        description = """Spectra Visual Editor is a tool for the treatment of spectrum files in order to identify lines, continuum regions and determine radial velocities among other options.
"""
        if "generate_spectrum" in dir(sve):
            description += """
The generation of synthetic spectrum is done thanks to:

SPECTRUM a Stellar Spectral Synthesis Program
(C) Richard O. Gray 1992 - 2010 Version 2.76e
"""

        licence = """Spectra Visual Editor is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Spectra Visual Editor is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Spectra Visual Editor.  If not, see <http://www.gnu.org/licenses/>."""


        info = wx.AboutDialogInfo()

        info.SetIcon(wx.Icon('images/SVE.png', wx.BITMAP_TYPE_PNG))
        info.SetName('Spectra Visual Editor')
        #info.SetVersion('2012.02.22')
        info.SetDescription(description)
        info.SetCopyright('(C) 2011 - 2012 Sergi Blanco Cuaresma')
        info.SetWebSite('http://www.marblestation.com')
        info.SetLicence(licence)
        info.AddDeveloper('Sergi Blanco Cuaresma')
        #info.AddDocWriter('Jan Bodnar')
        #info.AddArtist('The Tango crew')
        #info.AddTranslator('Jan Bodnar')

        wx.AboutBox(info)

    def flash_status_message(self, msg, flash_len_ms=3000, progress=True):
        if self.timeroff.IsRunning():
            self.timeroff.Stop()

        self.statusbar.SetStatusText(msg)
        if progress:
            self.gauge.SetValue(pos=100)

        self.Bind(
            wx.EVT_TIMER,
            self.on_flash_status_off,
            self.timeroff)
        self.timeroff.Start(flash_len_ms, oneShot=True)

    def status_message(self, msg):
        if self.timeroff.IsRunning():
            self.timeroff.Stop()
            # Progress bar to zero
            self.gauge.SetValue(pos=0)

        self.statusbar.SetStatusText(msg)

    def on_flash_status_off(self, event):
        self.statusbar.SetStatusText('')
        # Progress bar to zero
        self.gauge.SetValue(pos=0)


## Print usage
def usage():
    print "Usage:"
    print sys.argv[0], "[--continuum=file] [--lines=file] [--segments=file] [spectrum_file]"

## Interpret arguments
def get_arguments():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "cls", ["continuum=", "lines=", "segments="])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)

    continuum_file = None
    lines_file = None
    segments_file = None
    spectrum_file = None
    for o, a in opts:
        if o in ("-c", "--continuum"):
            continuum_file = a
            if not os.path.exists(continuum_file):
                print "Continuum file", args[0], "does not exists!"
                sys.exit(2)
        elif o in ("-l", "--lines"):
            lines_file = a
            if not os.path.exists(lines_file):
                print "Lines file", args[0], "does not exists!"
                sys.exit(2)
        elif o in ("-s", "--segments"):
            segments_file = a
            if not os.path.exists(segments_file):
                print "Segments file", args[0], "does not exists!"
                sys.exit(2)
        else:
            print "Argument", o, "not recognized!"
            usage()
            sys.exit(2)

    filenames = {}
    filenames['spectra'] = []
    filenames['continuum'] = continuum_file
    filenames['lines'] = lines_file
    filenames['segments'] = segments_file

    # Open spectra
    for arg in args:
        spectrum_file = arg
        if not os.path.exists(spectrum_file):
            print "Spectrum file", arg, "does not exists!"
            sys.exit(2)
        filenames['spectra'].append(spectrum_file)

    return filenames

#~ Example:
#~ python interactive.py --continuum=input/LUMBA/UVES_MRD_sun_cmask.txt --lines=input/LUMBA/UVES_MRD_sun_Fe-linelist.txt --segments=input/LUMBA/UVES_MRD_sun_segments.txt input/LUMBA/UVES_MRD_sun_official.s.gz

## Input files should be '\t' separated, comments with # and the first line is the header
if __name__ == '__main__':
    filenames = get_arguments()

    #### Read files

    spectra = []
    for path in filenames['spectra']:
        try:
            spectrum = sve.read_spectrum(path)
            #wfilter = (spectrum['waveobs'] >= 516.0) & (spectrum['waveobs'] <= 519.0)
            #spectrum = spectrum[wfilter]
        except Exception as e:
            print "Spectrum file", path, "has an incompatible format!"
            sys.exit(2)
        spectra.append(spectrum)

    if filenames['continuum'] != None:
        try:
            continuum = sve.read_continuum_regions(filenames['continuum'])
        except Exception as e:
            print "Continuum file", filenames['continuum'], "has an incompatible format!"
            sys.exit(2)

        ## Validations
        if np.any((continuum['wave_top'] - continuum['wave_base']) < 0):
            raise Exception("Continuum: wave_top cannot be smaller than wave_base")
    else:
        continuum = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])

    if filenames['lines'] != None:
        try:
            lines = sve.read_line_regions(filenames['lines'])
        except Exception as e:
            print "Lines file", filenames['lines'], "has an incompatible format!"
            sys.exit(2)


        ## Validations
        if np.any((lines['wave_top'] - lines['wave_base']) < 0):
            raise Exception("Lines: wave_top cannot be smaller than wave_base")

        if np.any((lines['wave_top'] - lines['wave_peak']) < 0):
            raise Exception("Lines: wave_top cannot be smaller than wave_peak")

        if np.any((lines['wave_peak'] - lines['wave_base']) < 0):
            raise Exception("Lines: wave_peak cannot be smaller than wave_base")
    else:
        lines = np.zeros((0,), dtype=[('wave_peak', '<f8'), ('wave_base', '<f8'), ('wave_top', '<f8')])

    if filenames['segments'] != None:
        try:
            segments = sve.read_segment_regions(filenames['segments'])
        except Exception as e:
            print "Segments file", filenames['segments'], "has an incompatible format!"
            sys.exit(2)

        ## Validations
        if np.any((segments['wave_top'] - segments['wave_base']) < 0):
            raise Exception("Segments: wave_top cannot be smaller than wave_base")
    else:
        segments = np.zeros((0,), dtype=[('wave_base', '<f8'), ('wave_top', '<f8')])



    regions = {}
    regions['continuum'] = continuum
    regions['lines'] = lines
    regions['segments'] = segments

    ## Force locale to English. If for instance Spanish is the default,
    ## wxPython makes that functions like atof(), called from external
    ## C libraries through Cython, behave interpreting comma as decimal separator
    ## and brokes the code that expects to have dots as decimal separators
    ## (i.e. SPECTRUM and its external files like linelists)
    locales = {
        u'en' : (wx.LANGUAGE_ENGLISH, u'en_US.UTF-8'),
        u'es' : (wx.LANGUAGE_SPANISH, u'es_ES.UTF-8'),
        u'fr' : (wx.LANGUAGE_FRENCH, u'fr_FR.UTF-8'),
        }
    os.environ['LANG'] = locales[u'en'][1]

    app = wx.App(False) # False to avoid additional windows for stdout/stderr
    app.frame = SpectraFrame(spectra, regions, filenames)
    app.frame.Show()
    app.MainLoop()


